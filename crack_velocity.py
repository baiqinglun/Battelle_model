"""
管道断裂力学计算公式
根据提供的四个公式实现相关计算
"""

import math
import matplotlib.pyplot as plt
import matplotlib

# 常量定义
COEFF_Vc = 0.670
COEFF_Pa = 0.382
COEFF_TERM = 3.81e7
COEFF_DpAp = 0.502

class CrackVelocityCalculator:
    """
    管道断裂力学计算类
    包含计算流动应力、裂纹止裂压力、裂纹速度和断裂阻力比的方法
    """

    def __init__(self, t = 6.3, D = 450.0, YS = 430, TS = 530, P = 8.0, Cv = 100.0, Ac = 100.0):
        self.t = t  # 壁厚 mm
        self.D = D  # 直径 mm
        self.YS = YS  # 屈服应力 MPa
        self.TS = TS  # 抗拉强度 MPa
        self.P = P  # 扩展压力 MPa
        self.Cv = Cv  # 夏比冲击能量 J
        self.Ac = Ac  # 韧带面积 mm²
        self.P_min = 1.0  # 最小压力 MPa
        self.P_max = 15.0  # 最大压力 MPa
        self.dP = 0.1
        
        self.results = {
            'pressure': [],
            'velocity': []
        }

    def calculate_flow_stress(self, YS, TS):
        """
        计算流动应力 σ_flow
        
        参数:
        YS: 屈服应力 (MPa)
        TS: 抗拉强度 (MPa)
        
        返回:
        σ_flow: 流动应力
        """
        return (YS + TS) / 2


    def calculate_arrest_pressure(self, t, D, sigma_flow, R):
        """
        计算裂纹止裂压力 Pa
        
        参数:
        t: 管道壁厚 (mm)
        D: 管道直径 (mm)
        sigma_flow: 流动应力
        R: 断裂阻力 (J/mm²)
        
        返回:
        Pa: 裂纹止裂压力 (MPa)
        """
        # 计算内部项
        inner_term = (-COEFF_TERM / math.sqrt(D * t)) * (R / (sigma_flow ** 2))
        
        # 计算指数和反余弦
        exp_term = math.exp(inner_term)
        arccos_term = math.acos(exp_term)  # math.acos返回弧度
        
        # 计算Pa
        Pa = COEFF_Pa * (t / D) * sigma_flow * arccos_term
        
        return Pa


    def calculate_crack_velocity(self, sigma_flow, R, P, Pa):
        """
        计算管道裂纹速度 Vc
        
        参数:
        sigma_flow: 流动应力
        R: 断裂阻力 (J/mm²)
        P: 裂纹尖端扩展压力 (MPa)
        Pa: 裂纹止裂压力 (MPa)
        
        返回:
        Vc: 管道裂纹速度 (m/s)
        """
        # 计算各项
        term1 = sigma_flow / math.sqrt(R)
        term2 = ((P / Pa) - 1) ** 0.393
        
        # 计算Vc
        Vc = COEFF_Vc * term1 * term2
        
        return Vc


    def calculate_fracture_resistance_ratio(self, t, Cv, Ac):
        """
        计算断裂阻力比 Dp/Ap
        
        参数:
        t: 管道壁厚 (mm)
        Cv: 夏比冲击吸收能量 (J)
        Ac: 预裂纹夏比试样韧带面积 (mm²)
        
        返回:
        DpAp: 断裂阻力比
        """
        # 计算各项
        term1 = t ** 0.5
        term2 = (Cv / Ac) ** 0.544
        
        # 计算Dp/Ap
        DpAp = COEFF_DpAp * term1 * term2
        
        return DpAp


    def calculate_R(self, t, Cv, Ac):
        """
        计算断裂阻力 R (直接计算版本)
        
        参数:
        t: 管道壁厚 (mm)
        Cv: 夏比冲击吸收能量 (J)
        Ac: 预裂纹夏比试样韧带面积 (mm²)
        
        返回:
        R: 断裂阻力 (J/mm²)
        """
        return self.calculate_fracture_resistance_ratio(t, Cv, Ac)


    def calculate_comprehensive(self):
        """
        综合计算所有参数
    
        参数:
        t: 管道壁厚 (mm)
        D: 管道直径 (mm)
        YS: 屈服应力 (MPa)
        TS: 抗拉强度 (MPa)
        P: 裂纹尖端扩展压力 (MPa)
        Cv: 夏比冲击吸收能量 (J)
        Ac: 预裂纹夏比试样韧带面积 (mm²)
        
        返回:
        dict: 包含所有计算结果的字典
        """
        # 1. 计算流动应力
        sigma_flow = self.calculate_flow_stress(self.YS, self.TS)
        
        # 2. 计算断裂阻力R
        R = self.calculate_fracture_resistance_ratio(self.t, self.Cv, self.Ac)
        
        # 3. 计算止裂压力Pa
        Pa = self.calculate_arrest_pressure(self.t, self.D, sigma_flow, R)
        
        # 4. 计算裂纹速度Vc
        Vc = self.calculate_crack_velocity(sigma_flow, R, self.P, Pa)
        
        # 5. 计算断裂阻力比
        DpAp = self.calculate_fracture_resistance_ratio(self.t, self.Cv, self.Ac)
        
        return {
            "sigma_flow": sigma_flow,
            "R": R,
            "Pa": Pa,
            "Vc": Vc,
            "Dp/Ap": DpAp
        }

    def calculate_decompression_velocity(self, P0, T0, gas_composition, gas_type="hydrogen"):
        """
        计算天然气管道的减压速度
        
        参数:
        P0: 初始压力 (MPa)
        T0: 初始温度 (K)
        gas_composition: 气体组分字典，如 {"CH4": 0.95, "C2H6": 0.03, ...}
        gas_type: 气体类型，可选 "natural_gas", "methane", "hydrogen" 等
        
        返回:
        a: 减压速度 (m/s)
        """
        # 不同气体的物理性质
        gas_properties = {
            "methane": {
                "gamma": 1.31,      # CH4的等熵指数
                "R": 518.3,         # CH4的气体常数 (J/(kg·K))
                "molar_mass": 16.04 # g/mol
            },
            "natural_gas": {
                "gamma": 1.28,      # 典型天然气
                "R": 518.3,         # 近似值
                "molar_mass": 18.0  # 近似平均分子量
            },
            "hydrogen": {
                "gamma": 1.41,
                "R": 4124.0,
                "molar_mass": 2.02
            }
        }
        
        # 获取气体性质
        if gas_type in gas_properties:
            props = gas_properties[gas_type]
        else:
            # 默认使用天然气性质
            props = gas_properties["natural_gas"]
        
        gamma = props["gamma"]
        R = props["R"]
        
        # 计算压缩因子Z（简化计算，实际中可能需要使用状态方程）
        # 对于高压天然气，Z通常小于1
        P_MPa = P0
        T_K = T0
        
        # 简化方法：使用理想气体近似，Z≈1
        # 或者使用经验公式计算Z
        if P_MPa < 5:  # 低压情况
            Z = 1.0
        elif P_MPa < 15:  # 中压情况
            Z = 0.95
        else:  # 高压情况
            Z = 0.9
        
        # 计算减压波速
        a = math.sqrt(gamma * R * T_K * Z)
        
        return a

    def calculate(self):
        """在给定参数下，绘制裂纹尖端扩展压力 P 与裂纹速度 Vc 的关系曲线。

        P 从 P_min 到 P_max（单位 MPa），步长为 step。
        """
        sigma_flow = self.calculate_flow_stress(self.YS, self.TS)
        R = self.calculate_fracture_resistance_ratio(self.t, self.Cv, self.Ac)
        Pa = self.calculate_arrest_pressure(self.t, self.D, sigma_flow, R)
        
        current_P = self.P_min
        while current_P <= self.P_max + 1e-8:
            self.results['pressure'].append(current_P)

            # 只有当 P ≥ Pa 时，(P/Pa - 1)**0.393 才在实数域有意义
            if current_P >= Pa:
                Vc = self.calculate_crack_velocity(sigma_flow, R, current_P, Pa)
            else:
                Vc = float("nan")  # 用 NaN 在曲线上产生间断，表示无物理意义区域

            self.results['velocity'].append(Vc)
            current_P += self.dP
            
        return self.results

        
    def print_initial_parameters(self):
        """打印初始参数信息"""
        print("管道断裂力学计算初始参数:")
        print(f"壁厚 t = {self.t} mm")
        print(f"直径 D = {self.D} mm")
        print(f"屈服应力 YS = {self.YS} MPa")
        print(f"抗拉强度 TS = {self.TS} MPa")
        print(f"扩展压力 P = {self.P} MPa")
        print(f"夏比冲击能量 Cv = {self.Cv} J")
        print(f"韧带面积 Ac = {self.Ac} mm²")
        print()

def main():
    CrackVelocityCalculator_instance = CrackVelocityCalculator()

    # 单独计算示例
    results = CrackVelocityCalculator_instance.calculate_comprehensive()
    for key, value in results.items():
        print(f"{key:10} = {value:.4f}")

    # 绘制初始压力 (1–15 MPa) 与裂纹速度的关系曲线
    CrackVelocityCalculator_instance.plot_velocity_vs_pressure(CrackVelocityCalculator_instance.t, CrackVelocityCalculator_instance.D, CrackVelocityCalculator_instance.YS, CrackVelocityCalculator_instance.TS, CrackVelocityCalculator_instance.Cv, CrackVelocityCalculator_instance.Ac,Tc = 300, P_min=1.0, P_max=15.0, step=1.0)
    
# 示例使用
# if __name__ == "__main__":
#     main()