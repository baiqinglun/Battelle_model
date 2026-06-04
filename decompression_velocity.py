"""
气体管道等熵流动模拟程序
基于以下三个公式：
1. u_i - u_{i-1} = (P_{i-1} - P_i) / (a_i * ρ_i)
2. c = √(r * P / ρ)  (声速公式，其中c就是a)
3. W_local = a - u
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
import math
# 方法A: 使用系统字体（推荐）
matplotlib.font_manager.fontManager.addfont("C:\\Windows\\Fonts\\simhei.ttf")
font_name = matplotlib.font_manager.FontProperties(fname="C:\\Windows\\Fonts\\simhei.ttf").get_name()
rcParams['font.sans-serif'] = [font_name]
rcParams['axes.unicode_minus'] = False

class IsentropicFlowSimulator:
    def __init__(self,gas_name = "hydrogen", P0=8e6, T0=293, u0=0.0, dt=1e-5, max_steps=10000, crack_v=200, A=0.4*0.05):
        """
        初始化等熵流动模拟器
        
        参数:
        gamma: 比热比 (空气约为1.4，天然气约为1.28-1.31)
        P0: 初始压力 (Pa)
        T0: 初始温度 (K)
        u0: 初始速度 (m/s)
        dt: 时间步长 (s)
        max_steps: 最大迭代步数
        """
        self.gas = {
            "hydrogen" : {
                "gamma": 1.41,
                "R": 4124.0,
                "molar_mass": 2.02
            },
            "methane" : {
                "gamma": 1.31,
                "R": 518.3,
                "molar_mass": 16.04
            },
        }
        
        self.gamma = self.gas[gas_name]["gamma"]
        self.R = self.gas[gas_name]["R"]
        self.molar = self.gas[gas_name]["molar_mass"]
        self.P0 = P0
        self.T0 = T0
        self.rho0 = self.P0 / self.R / self.T0
        self.u0 = u0
        self.dt = dt
        self.max_steps = max_steps
        self.crack_v = crack_v
        self.A = A
        self.t1 = 0.4 / crack_v
        self.t2 = 0.9 / crack_v
        self.L = 2.1
        self.D = 0.45
        self.V0 = math.pi * (self.D/2)**2 * self.L  # 初始体积 m³
        self.Pair = 0.0  # 停止压力 Pa
        self.R_ = 8.314  # J/(mol·K)
        
        # 存储模拟结果
        self.results = {
            'pressure': [],
            'velocity': [],
            'density': [],
            'sound_speed': [],
            'W_local': [],
            'time': [],
            'quantity': []
        }
    
    def calulate_leakage_area(self, time):
        """计算泄漏面积 A (m²)"""
        if time < self.t1:
            A = self.crack_v * time * 0.05
        elif time < self.t2:
            A = (0.4 + self.crack_v * (time - self.t1) * math.sqrt(2)) * (0.05 + self.crack_v * (time - self.t1) * math.sqrt(2))
        else:
            A = (0.4 + self.crack_v * (self.t2 - self.t1) * math.sqrt(2)) * (0.05 + self.crack_v * (self.t2 - self.t1) * math.sqrt(2))
        return A
    
    def calculate_sound_speed(self, P, rho):
        """计算声速: a = √(γ * P / ρ)"""
        return np.sqrt(self.gamma * P / rho)
    
    def calculate_W_local(self, a, u):
        """计算局部波速: W_local = a - u"""
        return a - u
    
    def calculate_Qv(self, time, P_current):
        temp = pow(2 / (self.gamma + 1), (self.gamma + 1)/(self.gamma - 1))
        v =  P_current * math.sqrt(self.gamma / (self.R * self.T0) * temp)
        print(f"计算得到的泄漏速度 v = {v:.3f} m/s")
        Qv = v * self.A
        print(f"计算得到的 Qv = {Qv:.6f} kg/s")
        return Qv

    def read_Qv_from_csv(self, time):
        """从CSV文件中读取Qv值"""
        import pandas as pd
        df = pd.read_csv('Q.csv')
        # 假设CSV文件有两列：time_ms和mass_flow_rate_kg_per_s
        # 将time_ms转换为秒
        df['time_s'] = df['time_ms'] / 1000.0
        # 找到最接近当前时间的行
        closest_row = df.iloc[(df['time_s'] - time).abs().argsort()[:1]]
        Qv = closest_row['mass_flow_rate_kg_per_s'].values[0]
        return Qv
    
    def update_state(self, P_prev, rho_prev, u_prev, m_prev, t):
        """
        使用公式1更新状态: u_i = u_{i-1} + (P_{i-1} - P_i)/(a_i * ρ_i)
        """
        Qv = self.calculate_Qv(t, P_prev)
        Qm = Qv * self.dt * rho_prev
        m_current = m_prev - Qm  # 质量守恒
        P_current = P_prev * m_current / m_prev
        rho_current = rho_prev * (P_current / P_prev) ** (1/self.gamma)
        
        # 计算新的声速
        a_new = self.calculate_sound_speed(P_current, rho_current)
        
        # 使用公式1计算速度变化
        delta_u = (P_prev - P_current) / (a_new * rho_current)
        u_current = u_prev + delta_u
        
        return P_current, rho_current, u_current, a_new, m_current
    
    def simulate(self, target_W=0.0, tolerance=1e-6):
        """
        运行模拟，直到W_local接近0
        参数:
        target_W: 目标W_local值
        tolerance: 允许的误差
        """
        # 初始化状态
        P = self.P0
        rho = self.rho0
        u = self.u0
        a = self.calculate_sound_speed(P, rho)
        W = self.calculate_W_local(a, u)
        m = rho * self.V0  # 初始质量
        
        t = 0.0
        step = 0
        
        print("开始等熵流动模拟...")
        print(f"初始条件: P={P:.1f} Pa, ρ={rho:.3f} kg/m³, u={u:.3f} m/s, a={a:.2f} m/s, W={W:.2f} m/s")
        print("-" * 80)
        
        while step < self.max_steps:
            # 保存当前状态
            self.results['pressure'].append(P)
            self.results['velocity'].append(u)
            self.results['density'].append(rho)
            self.results['sound_speed'].append(a)
            self.results['W_local'].append(W)
            self.results['time'].append(t)
            self.results['quantity'].append(rho * self.V0)  # 假设体积V0不变
            
            # 检查是否达到目标
            if np.abs(W - target_W) < tolerance:
                print(f"\n在t={t:.6f}s时达到目标: W_local = {W:.6f} m/s ≈ {target_W} m/s")
                print(f"最终状态: P={P:.1f} Pa, u={u:.3f} m/s")
                break
            if W < 0:
                print(f"\n在t={t:.6f}s时W_local变为负值: W_local = {W:.6f} m/s")
                print(f"最终状态: P={P:.1f} Pa, u={u:.3f} m/s")
                print("+" * 80)
                break
            
            # 更新状态
            P, rho, u, a, m = self.update_state(P, rho, u, m, t)
            W = self.calculate_W_local(a, u)
            
            t += self.dt
            step += 1
            
            # 每1000步显示进度
            if step % 1000 == 0:
                print(f"步数: {step:5d}, 时间: {t:.6f}s, W={W:.3f} m/s, P={P:.1f} Pa, u={u:.3f} m/s")
        
        if step == self.max_steps:
            print(f"\n达到最大步数 {self.max_steps}，W_local = {W:.6f} m/s")
        
        return self.results
    
    def plot_results(self):
        """绘制模拟结果"""
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # 1. 压力-速度曲线
        ax1 = axes[0, 0]
        ax1.plot(self.results['velocity'], self.results['pressure'], 'b-', linewidth=2)
        ax1.set_xlabel('流速 u (m/s)')
        ax1.set_ylabel('压力 P (Pa)')
        ax1.set_title('压力-速度关系曲线')
        ax1.grid(True, alpha=0.3)
        
        # 2. W_local随时间变化
        ax2 = axes[0, 1]
        ax2.plot(self.results['time'], self.results['W_local'], 'r-', linewidth=2)
        ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='W=0')
        ax2.set_xlabel('时间 t (s)')
        ax2.set_ylabel('局部波速 W_local (m/s)')
        ax2.set_title('局部波速随时间变化')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. 压力随时间变化
        ax3 = axes[0, 2]
        ax3.plot(self.results['time'], self.results['pressure'], 'g-', linewidth=2)
        ax3.set_xlabel('时间 t (s)')
        ax3.set_ylabel('压力 P (Pa)')
        ax3.set_title('压力随时间变化')
        ax3.grid(True, alpha=0.3)
        
        # 4. 速度随时间变化
        ax4 = axes[1, 0]
        ax4.plot(self.results['time'], self.results['velocity'], 'm-', linewidth=2)
        ax4.set_xlabel('时间 t (s)')
        ax4.set_ylabel('流速 u (m/s)')
        ax4.set_title('速度随时间变化')
        ax4.grid(True, alpha=0.3)
        
        # 5. 声速随时间变化
        ax5 = axes[1, 1]
        ax5.plot(self.results['time'], self.results['sound_speed'], 'c-', linewidth=2)
        ax5.set_xlabel('时间 t (s)')
        ax5.set_ylabel('声速 a (m/s)')
        ax5.set_title('声速随时间变化')
        ax5.grid(True, alpha=0.3)
        
        # 6. 密度随时间变化
        ax6 = axes[1, 2]
        ax6.plot(self.results['time'], self.results['density'], 'y-', linewidth=2)
        ax6.set_xlabel('时间 t (s)')
        ax6.set_ylabel('密度 ρ (kg/m³)')
        ax6.set_title('密度随时间变化')
        ax6.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # 单独绘制相图：W_local vs (a-u)
        fig2, ax = plt.subplots(figsize=(10, 6))
        
        # 计算a-u
        a_minus_u = [a - u for a, u in zip(self.results['sound_speed'], self.results['velocity'])]
        
        ax.plot(self.results['time'], a_minus_u, 'b-', linewidth=2, label='a - u')
        ax.plot(self.results['time'], self.results['W_local'], 'r--', linewidth=2, label='W_local')
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='零线')
        
        ax.set_xlabel('时间 t (s)')
        ax.set_ylabel('速度差 (m/s)')
        ax.set_title('声速与流速差随时间变化')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # 单独绘制相图：W_local vs (a-u)
        fig3, ax = plt.subplots(figsize=(10, 6))
        
        # 计算a-u
        a_minus_u = [a - u for a, u in zip(self.results['sound_speed'], self.results['velocity'])]
        ax.plot(a_minus_u, self.results['pressure'], 'b-', linewidth=2, label='a - u')
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.5, label='零线')
        
        ax.set_xlabel('减压速度 u (m/s)')
        ax.set_ylabel('压力 P (Pa)')
        ax.set_title('压力-速度关系曲线')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        return fig, fig2
    
    def print_summary(self):
        """打印模拟结果摘要"""
        print("\n" + "="*60)
        print("模拟结果摘要")
        print("="*60)
        print(f"总模拟时间: {self.results['time'][-1]:.6f} s")
        print(f"总步数: {len(self.results['time'])}")
        print(f"初始压力: {self.results['pressure'][0]:.1f} Pa")
        print(f"最终压力: {self.results['pressure'][-1]:.1f} Pa")
        print(f"压力变化: {self.results['pressure'][0] - self.results['pressure'][-1]:.1f} Pa")
        print(f"初始速度: {self.results['velocity'][0]:.3f} m/s")
        print(f"最终速度: {self.results['velocity'][-1]:.3f} m/s")
        print(f"速度变化: {self.results['velocity'][-1] - self.results['velocity'][0]:.3f} m/s")
        print(f"初始W_local: {self.results['W_local'][0]:.3f} m/s")
        print(f"最终W_local: {self.results['W_local'][-1]:.6f} m/s")
        print("="*60)


def main():
    """主函数：运行模拟并绘制结果"""
    # 创建模拟器实例
    # 参数说明：
    simulator = IsentropicFlowSimulator(
        gas_name="hydrogen",
        P0=8e6,      # 800 kPa
        T0=293,      # 293 K
        u0=0.0,      # 初始静止
        dt=1e-6,     # 时间步长
        max_steps=1e5,
        crack_v=200
    )
    
    # 运行模拟，直到W_local接近0
    results = simulator.simulate(target_W=0.0, tolerance=1e-3)
    
    # 打印摘要
    simulator.print_summary()
    
    # 绘制图表
    simulator.plot_results()
    
    # 额外分析：找出W_local=0附近的点
    W_array = np.array(results['W_local'])
    idx_near_zero = np.argmin(np.abs(W_array))
    
    print(f"\nW_local最接近0的点:")
    print(f"  时间: {results['time'][idx_near_zero]:.6f} s")
    print(f"  W_local: {results['W_local'][idx_near_zero]:.6f} m/s")
    print(f"  压力: {results['pressure'][idx_near_zero]:.1f} Pa")
    print(f"  速度: {results['velocity'][idx_near_zero]:.3f} m/s")
    print(f"  声速: {results['sound_speed'][idx_near_zero]:.3f} m/s")
    print(f"  a - u: {results['sound_speed'][idx_near_zero] - results['velocity'][idx_near_zero]:.6f} m/s")


if __name__ == "__main__":
    main()