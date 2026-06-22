from decompression_velocity import IsentropicFlowSimulator
from crack_velocity import CrackVelocityCalculator
from tools import draw

hydrogen_simulator = IsentropicFlowSimulator(
        gas_name="hydrogen",
        P0=8.25e6,      # 8MPa
        T0=300,    # 293 K
        u0=0.0,      # 初始静止
        dt=1e-4,     # 时间步长
        max_steps=1e4,
        crack_v = 200,
        A = 0.1516
    )
hydrogen_results = hydrogen_simulator.simulate(target_W=0.0, tolerance=1e-3)

CrackVelocityCalculator_instance1 = CrackVelocityCalculator(
    t=6.3,        # 管道壁厚 (mm)
    D=463.0,       # 管道直径 (mm)
    YS=430.0,      # 屈服应力 (MPa)
    TS=530.0,      # 抗拉强度 (MPa)
    Pmin=1,         # 最小裂纹尖端扩展压力 (MPa)
    Pmax=8.25,        # 最大裂纹尖端扩展压力 (MPa)
    Cv=50,      # 夏比冲击吸收能量 (J)
    Ac=80.0      # 预裂纹夏比试样韧带面积 (mm²)
)
crack_results = CrackVelocityCalculator_instance1.calculate()

# 绘图
hydrogen_a_minus_u = [a - u for a, u in zip(hydrogen_results['sound_speed'], hydrogen_results['velocity'])]
hydrogen_scaled_gas_pressure = [p / 1e6 for p in hydrogen_results['pressure']]
scaled_crack_pressure = [p for p in crack_results['pressure']]
draw([{"x":hydrogen_a_minus_u, "y":hydrogen_scaled_gas_pressure, "style":'b-', "label":'hydrogen'},
    #   {"x":methane_a_minus_u, "y":methane_scaled_gas_pressure, "style":'g-', "label":'methane'},
      {"x":crack_results['velocity'], "y":scaled_crack_pressure, "style":'r--', "label":'crack'}],
     xlabel='减压速度 u (m/s)',
     ylabel='压力 P (MPa)',
     title='压力-速度关系曲线',
     xlim=(0, 1500),
     ylim=(0, 12))

# 绘制时间-流量和时间-压降曲线
# hydrogen_simulator.plot_flow_and_pressure_drop()
