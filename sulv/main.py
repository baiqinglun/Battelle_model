import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, PchipInterpolator
from scipy.integrate import cumulative_trapezoid
import warnings
import math
import pandas as pd
from pathlib import Path
from scipy.signal import savgol_filter
warnings.filterwarnings('ignore', category=DeprecationWarning)

plt.rcParams['font.family'] = 'SimHei'  # 替换为你选择的字体

def cumtrapz(y, x=None, dx=1.0, initial=None):
    """
    cumtrapz的兼容性函数
    在SciPy >= 1.6.0中，cumtrapz已被cumulative_trapezoid取代
    """
    try:
        # 尝试使用较新版本的cumulative_trapezoid
        from scipy.integrate import cumulative_trapezoid
        result = cumulative_trapezoid(y, x=x, dx=dx, initial=initial)
    except (ImportError, AttributeError):
        # 如果不存在，使用numpy实现
        if x is None:
            x = np.arange(len(y)) * dx
        
        if len(y) < 2:
            return np.array([0])
        
        # 使用梯形法则手动计算累积积分
        result = np.zeros(len(y))
        for i in range(1, len(y)):
            dx_i = x[i] - x[i-1]
            result[i] = result[i-1] + 0.5 * (y[i] + y[i-1]) * dx_i
        
        if initial is not None:
            result += initial
    
    return result

class PipelineLeakageAnalyzer:
    """管道泄漏分析器"""
    
    def __init__(self, gas_type='air', isothermal=True):
        """
        初始化分析器
        
        参数:
        ----------
        gas_type : str
            气体类型 ('air', 'natural_gas', 'methane' 等)
        isothermal : bool
            是否为等温过程
        """
        # 设置气体参数
        self.gas_properties = {
            'air': {
                'gamma': 1.4,           # 比热比
                'R': 287.058,          # 气体常数 [J/(kg·K)]
                'M': 0.02897,          # 摩尔质量 [kg/mol]
                'cp': 1005,            # 定压比热 [J/(kg·K)]
                'cv': 718              # 定容比热 [J/(kg·K)]
            },
            'hydrogen': {
                'gamma': 1.41,         # 比热比
                'R': 4157.0,           # 气体常数 [J/(kg·K)]
                'M': 0.002,            # 摩尔质量 [kg/mol]
                'cp': 14300,           # 定压比热 [J/(kg·K)]
                'cv': 10100            # 定容比热 [J/(kg·K)]
            },
            'natural_gas': {
                'gamma': 1.31,
                'R': 518.3,
                'M': 0.01604,
                'cp': 2200,
                'cv': 1680
            },
            'methane': {
                'gamma': 1.32,
                'R': 518.3,
                'M': 0.01604,
                'cp': 2250,
                'cv': 1700
            }
        }
        
        self.gas = self.gas_properties.get(gas_type, self.gas_properties['air'])
        self.isothermal = isothermal
        self.gamma = self.gas['gamma']
        self.R = self.gas['R']

    def smooth_pressure(self, pressure_pa, window_length=51, polyorder=2):
        """对压力序列进行平滑处理"""
        pressure_pa = np.asarray(pressure_pa, dtype=float)

        if pressure_pa.size < 3:
            return pressure_pa.copy()

        window_length = min(window_length, pressure_pa.size if pressure_pa.size % 2 == 1 else pressure_pa.size - 1)
        if window_length < 3:
            return pressure_pa.copy()

        if window_length % 2 == 0:
            window_length -= 1

        polyorder = min(polyorder, window_length - 1)
        if polyorder < 1:
            return pressure_pa.copy()

        try:
            smoothed = savgol_filter(pressure_pa, window_length=window_length, polyorder=polyorder, mode='interp')
        except Exception:
            kernel = np.ones(window_length) / window_length
            smoothed = np.convolve(pressure_pa, kernel, mode='same')

        # 压力应当随时间非增，避免平滑带来的局部回升
        return np.minimum.accumulate(smoothed)

    def smooth_derivative(self, pressure_pa, time, window_length=51, polyorder=2):
        """直接求平滑导数，避免两次平滑"""
        pressure_pa = np.asarray(pressure_pa, dtype=float)
        time = np.asarray(time, dtype=float)

        if pressure_pa.size < 3:
            return np.zeros_like(pressure_pa)

        window_length = min(window_length, pressure_pa.size if pressure_pa.size % 2 == 1 else pressure_pa.size - 1)
        if window_length < 3:
            return np.gradient(pressure_pa, time)

        if window_length % 2 == 0:
            window_length -= 1

        polyorder = min(polyorder, window_length - 1)
        if polyorder < 1:
            return np.gradient(pressure_pa, time)

        try:
            # 使用 Savitzky-Golay 的 deriv 参数直接求一阶导数，避免两次平滑
            dt = np.mean(np.diff(time))
            dp_dt = savgol_filter(pressure_pa, window_length=window_length, polyorder=polyorder, deriv=1, delta=dt, mode='interp')
            return dp_dt
        except Exception:
            # 降级方案：对平滑后的压力求导
            smoothed = self.smooth_pressure(pressure_pa, window_length=window_length, polyorder=polyorder)
            return np.gradient(smoothed, time)
        
    def calculate_leakage_rate(self, time_ms, pressure_mpa, A_leak, V_pipe, 
                               P0=None, T0=300.0, P_back=0.101325,
                               smooth_window=51, smooth_polyorder=2):
        """
        计算泄漏速率
        
        参数:
        ----------
        time_ms : array-like
            时间数组 [ms]
        pressure_mpa : array-like
            压力数组 [MPa]
        A_leak : float
            泄漏面积 [m²]
        V_pipe : float
            管道容积 [m³]
        P0 : float, optional
            初始压力 [MPa], 如果为None则使用第一个压力值
        T0 : float
            初始温度 [K]
        P_back : float
            背压 [MPa], 默认大气压
        
        返回:
        ----------
        dict: 包含各种计算结果
        """
        # 单位转换
        time = np.array(time_ms) * 1e-3  # ms -> s
        pressure = np.array(pressure_mpa)  # MPa
        
        # 转换压力为Pa
        pressure_raw_pa = pressure * 1e6
        P_back_pa = P_back * 1e6

        pressure_pa = self.smooth_pressure(
            pressure_raw_pa,
            window_length=smooth_window,
            polyorder=smooth_polyorder
        )
        
        if P0 is None:
            P0 = pressure_pa[0]
        else:
            P0 = P0 * 1e6
        
        # 1. 压力对时间的导数 (dp/dt) - 使用平滑导数
        dp_dt = self.smooth_derivative(pressure_pa, time, window_length=smooth_window, polyorder=smooth_polyorder)
        
        # 2. 插值函数用于平滑
        f_pressure = PchipInterpolator(time, pressure_pa, extrapolate=True)
        time_fine = np.linspace(time.min(), time.max(), 1000)
        pressure_fine = f_pressure(time_fine)
        dp_dt_fine = np.gradient(pressure_fine, time_fine)
        
        # 3. 使用质量守恒法计算质量泄漏率
        # 对于理想气体: ρ = P/(R*T)
        if self.isothermal:
            # 等温过程: d(ρ)/dt = (1/(R*T)) * dP/dt
            rho = pressure_pa / (self.R * T0)  # 密度 [kg/m³]
            drho_dt = dp_dt / (self.R * T0)
            
            # 质量泄漏率: dm/dt = -V * dρ/dt
            mass_flow_rate = -V_pipe * drho_dt  # [kg/s]
            mass_flow_rate = np.maximum(mass_flow_rate, 0)
            
            # 体积流量: Q = dm/dt / ρ
            volumetric_flow_rate = mass_flow_rate / rho  # [m³/s]
            
        else:
            # 绝热过程: 使用更复杂的公式
            # 假设为理想气体，绝热膨胀
            mass_flow_rate = np.zeros_like(time)
            for i in range(len(time)):
                P = pressure_pa[i]
                if P > P_back_pa:
                    # 判断是否为临界流
                    P_ratio = P_back_pa / P
                    critical_pressure_ratio = (2/(self.gamma+1))**(self.gamma/(self.gamma-1))
                    
                    if P_ratio <= critical_pressure_ratio:
                        # 临界流（壅塞流）
                        mass_flow_rate[i] = A_leak * P * np.sqrt(
                            self.gamma/(self.R*T0) * (2/(self.gamma+1))**((self.gamma+1)/(self.gamma-1))
                        )
                    else:
                        # 亚临界流
                        mass_flow_rate[i] = A_leak * P * np.sqrt(
                            2*self.gamma/((self.gamma-1)*self.R*T0) * 
                            (P_ratio**(2/self.gamma) - P_ratio**((self.gamma+1)/self.gamma))
                        )
                else:
                    mass_flow_rate[i] = 0
            
            volumetric_flow_rate = mass_flow_rate * self.R * T0 / pressure_pa
        
        # 4. 使用孔口流动公式验证
        # 等熵孔口流动公式
        orifice_flow_rate = np.zeros_like(pressure_pa)
        for i, P in enumerate(pressure_pa):
            if P > P_back_pa:
                beta = P_back_pa / P
                critical_beta = (2/(self.gamma+1))**(self.gamma/(self.gamma-1))
                
                if beta <= critical_beta:
                    # 壅塞流
                    orifice_flow_rate[i] = A_leak * P * np.sqrt(
                        self.gamma/(self.R*T0) * (2/(self.gamma+1))**((self.gamma+1)/(self.gamma-1))
                    )
                else:
                    # 亚壅塞流
                    orifice_flow_rate[i] = A_leak * P * np.sqrt(
                        (2*self.gamma/((self.gamma-1)*self.R*T0)) * 
                        (beta**(2/self.gamma) - beta**((self.gamma+1)/self.gamma))
                    )
            else:
                orifice_flow_rate[i] = 0
        orifice_flow_rate = np.maximum(orifice_flow_rate, 0)

        if A_leak > 0:
            leak_velocity = volumetric_flow_rate / A_leak
        else:
            leak_velocity = np.zeros_like(volumetric_flow_rate)
        
        # 5. 计算累积泄漏量
        cumulative_leak = cumtrapz(mass_flow_rate, time, initial=0)
        
        # 6. 计算特征时间
        # 找到压力下降到初始压力的1/e的时间
        P_target = P0 / np.e
        idx = np.where(pressure_pa <= P_target)[0]
        if len(idx) > 0:
            tau = time[idx[0]]
        else:
            tau = time[-1]
        
        return {
            'time': time,
            'time_fine': time_fine,
            'pressure_raw': pressure_raw_pa,
            'pressure': pressure_pa,
            'pressure_fine': pressure_fine,
            'dp_dt': dp_dt,
            'dp_dt_fine': dp_dt_fine,
            'mass_flow_rate': mass_flow_rate,
            'volumetric_flow_rate': volumetric_flow_rate,
            'leak_velocity': leak_velocity,
            'orifice_flow_rate': orifice_flow_rate,
            'cumulative_leak': cumulative_leak,
            'density': rho if self.isothermal else None,
            'time_constant': tau,
            'parameters': {
                'A_leak': A_leak,
                'V_pipe': V_pipe,
                'P0': P0,
                'T0': T0,
                'P_back': P_back_pa,
                'gas_properties': self.gas
            }
        }
    
    def plot_results(self, results, figsize=(15, 10)):
        """绘制结果"""
        fig, axes = plt.subplots(3, 2, figsize=figsize)
        
        # 1. 压力-时间曲线
        ax1 = axes[0, 0]
        ax1.scatter(results['time']*1000, results['pressure_raw']/1e6, s=16, alpha=0.5, label='原始数据')
        ax1.plot(results['time']*1000, results['pressure']/1e6, 'b-', linewidth=2, label='平滑曲线')
        ax1.plot(results['time_fine']*1000, results['pressure_fine']/1e6, 'r--', 
            alpha=0.7, label='插值曲线')
        ax1.set_xlabel('时间 (ms)')
        ax1.set_ylabel('压力 (MPa)')
        ax1.set_title('压力衰减曲线')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 2. 压力变化率
        ax2 = axes[0, 1]
        ax2.plot(results['time']*1000, results['dp_dt']/1e6, 'g-', linewidth=2)
        ax2.set_xlabel('时间 (ms)')
        ax2.set_ylabel('dp/dt (MPa/s)')
        ax2.set_title('压力变化率')
        ax2.grid(True, alpha=0.3)
        
        # 3. 质量泄漏率
        ax3 = axes[1, 0]
        ax3.plot(results['time']*1000, results['mass_flow_rate'], 'r-', 
                linewidth=2, label='质量守恒法')
        ax3.plot(results['time']*1000, results['orifice_flow_rate'], 'b--', 
                linewidth=2, alpha=0.7, label='孔口公式')
        ax3.set_xlabel('时间 (ms)')
        ax3.set_ylabel('质量泄漏率 (kg/s)')
        ax3.set_title('泄漏速率对比')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        # 4. 体积流量
        ax4 = axes[1, 1]
        ax4.plot(results['time']*1000, results['volumetric_flow_rate'], 'purple', linewidth=2)
        ax4.set_xlabel('时间 (ms)')
        ax4.set_ylabel('体积流量 (m³/s)')
        ax4.set_title('体积泄漏流量')
        ax4.grid(True, alpha=0.3)
        
        # 5. 累积泄漏量
        ax5 = axes[2, 0]
        ax5.plot(results['time']*1000, results['cumulative_leak'], 'brown', linewidth=2)
        ax5.set_xlabel('时间 (ms)')
        ax5.set_ylabel('累积泄漏量 (kg)')
        ax5.set_title('累积泄漏质量')
        ax5.grid(True, alpha=0.3)
        
        # 6. 相图: 泄漏率 vs 压力
        ax6 = axes[2, 1]
        scatter = ax6.scatter(results['pressure']/1e6, results['mass_flow_rate'], 
                            c=results['time'], cmap='viridis', s=20)
        ax6.set_xlabel('压力 (MPa)')
        ax6.set_ylabel('质量泄漏率 (kg/s)')
        ax6.set_title('泄漏率-压力相图')
        ax6.grid(True, alpha=0.3)
        plt.colorbar(scatter, ax=ax6, label='时间 (s)')
        
        plt.tight_layout()
        plt.show()
        
        return fig

    def plot_leakage_velocity(self, results, figsize=(10, 5)):
        """绘制泄漏速度曲线"""
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(results['time'] * 1000, results['leak_velocity'], color='darkorange', linewidth=2)
        ax.set_xlabel('时间 (ms)')
        ax.set_ylabel('泄漏速度 (m/s)')
        ax.set_title('泄漏速度曲线')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

        return fig

    def export_results(self, results, output_file='leakage_results.csv'):
        """导出计算结果到CSV文件"""
        export_data = {
            'time_ms': results['time'] * 1000,
            'pressure_mpa': results['pressure'] / 1e6,
            'pressure_raw_mpa': results['pressure_raw'] / 1e6,
            'dp_dt_mpa_per_s': results['dp_dt'] / 1e6,
            'mass_flow_rate_kg_per_s': results['mass_flow_rate'],
            'volumetric_flow_rate_m3_per_s': results['volumetric_flow_rate'],
            'leak_velocity_m_per_s': results['leak_velocity'],
            'orifice_flow_rate_kg_per_s': results['orifice_flow_rate'],
            'cumulative_leak_kg': results['cumulative_leak']
        }
        
        df = pd.DataFrame(export_data)
        df.to_csv(output_file, index=False, float_format='%.6f')
        print(f"\n结果已导出到: {output_file}")
    
    def print_summary(self, results):
        """打印结果摘要"""
        params = results['parameters']
        
        print("="*60)
        print("管道泄漏分析结果摘要")
        print("="*60)
        print(f"泄漏面积: {params['A_leak']:.6f} m²")
        print(f"管道容积: {params['V_pipe']:.3f} m³")
        print(f"初始压力: {params['P0']/1e6:.3f} MPa")
        print(f"初始温度: {params['T0']:.1f} K")
        print(f"背压: {params['P_back']/1e6:.3f} MPa")
        print(f"气体类型: R={self.gas['R']:.1f} J/(kg·K), γ={self.gamma:.3f}")
        print("-"*60)
        
        # 计算统计量
        mfr = results['mass_flow_rate']
        peak_flow = np.max(mfr)
        avg_flow = np.mean(mfr)
        total_leak = results['cumulative_leak'][-1]
        
        print(f"峰值泄漏率: {peak_flow:.4f} kg/s")
        print(f"平均泄漏率: {avg_flow:.4f} kg/s")
        print(f"总泄漏质量: {total_leak:.4f} kg")
        print(f"特征时间常数: {results['time_constant']*1000:.2f} ms")
        print("="*60)


# 示例使用代码
if __name__ == "__main__":
    # 读取csv文件
    data_path = Path(__file__).resolve().parent / "6_data.csv"
    data = pd.read_csv(data_path)
    time_ms = data['time_ms'].values
    pressure_mpa = data['pressure_mpa'].values

    # 已知参数
    A_leak = 0.0091  # 泄漏面积 [m²] (来自图片)
    
    V_pipe = math.pi * (462.6 / 1000 / 2) ** 2 * 1  # 管道容积 [m³]，按毫米换算为米
    
    # 创建分析器实例
    analyzer = PipelineLeakageAnalyzer(gas_type='hydrogen', isothermal=True)
    
    # 计算泄漏速率
    results = analyzer.calculate_leakage_rate(
        time_ms=time_ms,
        pressure_mpa=pressure_mpa,
        A_leak=A_leak,
        V_pipe=V_pipe,
        P0=None,  # 初始压力 [MPa]，自动使用数据起点
        T0=300.0,  # 温度 [K]
        P_back=0.101325,  # 大气压 [MPa]
        smooth_window=51,
        smooth_polyorder=2
    )
    
    # 打印结果摘要
    analyzer.print_summary(results)
    
    # 绘制图形
    analyzer.plot_results(results)
    analyzer.plot_leakage_velocity(results)
    
    # 输出详细数据
    print("\n详细数据表:")
    print("时间(ms)  压力(MPa)  质量泄漏率(kg/s)  泄漏速度(m/s)  累积泄漏(kg)")
    print("-"*75)
    for i in range(len(results['time'])):
        t = results['time'][i] * 1000
        p = results['pressure'][i] / 1e6
        mfr = results['mass_flow_rate'][i]
        vel = results['leak_velocity'][i]
        cum = results['cumulative_leak'][i]
        print(f"{t:6.2f}   {p:8.3f}   {mfr:15.6f}   {vel:12.6f}   {cum:12.6f}")
    
    # 导出结果到CSV
    output_file = data_path.parent / 'leakage_results.csv'
    analyzer.export_results(results, str(output_file))


# 更高级的分析：使用曲线拟合优化参数
def optimize_leakage_analysis(time_ms, pressure_mpa, A_leak, V_pipe_guess=0.1):
    """
    通过拟合优化管道容积参数
    
    参数:
    ----------
    time_ms, pressure_mpa : 实验数据
    A_leak : 泄漏面积
    V_pipe_guess : 管道容积的初始猜测值
    
    返回:
    ----------
    优化的容积和拟合结果
    """
    from scipy.optimize import curve_fit
    
    def pressure_model(t, V, P0, tau):
        """指数衰减模型: P = P0 * exp(-t/tau)"""
        return P0 * np.exp(-t / tau)
    
    # 准备数据
    time_s = np.array(time_ms) * 1e-3
    pressure_pa = np.array(pressure_mpa) * 1e6
    
    # 移除零压力点
    mask = pressure_pa > 0.1e6
    time_fit = time_s[mask]
    pressure_fit = pressure_pa[mask]
    
    if len(time_fit) < 3:
        print("数据点不足，无法进行拟合")
        return V_pipe_guess, None
    
    # 初始猜测
    p0_guess = pressure_fit[0]
    tau_guess = time_fit[-1] / 3  # 初始时间常数猜测
    
    try:
        # 拟合曲线
        popt, pcov = curve_fit(pressure_model, time_fit, pressure_fit, 
                              p0=[p0_guess, tau_guess], maxfev=5000)
        
        P0_fit, tau_fit = popt
        
        # 从时间常数估算容积
        # 对于等温过程: τ = V / (A_leak * sqrt(RT))
        analyzer = PipelineLeakageAnalyzer()
        R = analyzer.R
        T0 = 300.0
        V_fit = A_leak * np.sqrt(R * T0) * tau_fit
        
        print("\n" + "="*60)
        print("曲线拟合优化结果")
        print("="*60)
        print(f"拟合初始压力: {P0_fit/1e6:.3f} MPa")
        print(f"拟合时间常数: {tau_fit*1000:.2f} ms")
        print(f"估算管道容积: {V_fit:.6f} m³")
        print("="*60)
        
        return V_fit, {'P0_fit': P0_fit, 'tau_fit': tau_fit, 'V_fit': V_fit}
        
    except Exception as e:
        print(f"拟合失败: {e}")
        return V_pipe_guess, None


# 使用优化函数
print("\n\n进行参数优化...")
optimized_V, fit_results = optimize_leakage_analysis(
    time_ms=time_ms,
    pressure_mpa=pressure_mpa,
    A_leak=A_leak,
    V_pipe_guess=0.1
)

# 如果优化成功，使用优化后的容积重新计算
if fit_results is not None and optimized_V > 0:
    print(f"\n使用优化后的容积 {optimized_V:.6f} m³ 重新计算...")
    analyzer2 = PipelineLeakageAnalyzer()
    results2 = analyzer2.calculate_leakage_rate(
        time_ms=time_ms,
        pressure_mpa=pressure_mpa,
        A_leak=A_leak,
        V_pipe=optimized_V,
        P0=fit_results['P0_fit']/1e6,
        T0=300.0
    )
    analyzer2.print_summary(results2)