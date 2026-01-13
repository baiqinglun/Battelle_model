from decompression_velocity import IsentropicFlowSimulator
from crack_velocity import CrackVelocityCalculator
from tools import draw

simulator = IsentropicFlowSimulator(
        gas_name="hydrogen",
        P0=8e6,      # 8MPa
        T0=293,    # 293 K
        u0=0.0,      # 初始静止
        dt=1e-5,     # 时间步长
        max_steps=1e5
    )
gas_results = simulator.simulate(target_W=0.0, tolerance=1e-3)

CrackVelocityCalculator_instance = CrackVelocityCalculator(
    t=6,        # 管道壁厚 (mm)
    D=450.0,       # 管道直径 (mm)
    YS=400.0,      # 屈服应力 (MPa)
    TS=500.0,      # 抗拉强度 (MPa)
    Pmin=1,         # 最小裂纹尖端扩展压力 (MPa)
    Pmax=12,        # 最大裂纹尖端扩展压力 (MPa)
    Cv=50.0,      # 夏比冲击吸收能量 (J)
    Ac=80.0      # 预裂纹夏比试样韧带面积 (mm²)
)
crack_results = CrackVelocityCalculator_instance.calculate()

# 绘图
a_minus_u = [a - u for a, u in zip(gas_results['sound_speed'], gas_results['velocity'])]
scaled_gas_pressure = [p / 1e6 for p in gas_results['pressure']]
scaled_crack_pressure = [p for p in crack_results['pressure']]
draw([{"x":a_minus_u, "y":scaled_gas_pressure, "style":'b-', "label":'gas'},
      {"x":crack_results['velocity'], "y":scaled_crack_pressure, "style":'r--', "label":'crack'}],
     xlabel='减压速度 u (m/s)',
     ylabel='压力 P (MPa)',
     title='压力-速度关系曲线',
     xlim=(0, 1500),
     ylim=(0, 12))
