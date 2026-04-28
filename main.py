from decompression_velocity import IsentropicFlowSimulator
from crack_velocity import CrackVelocityCalculator
from tools import draw

hydrogen_simulator = IsentropicFlowSimulator(
        gas_name="hydrogen",
        P0=8.24e6,      # 8MPa
        T0=300,    # 293 K
        u0=0.0,      # 初始静止
        dt=1e-4,     # 时间步长
        max_steps=1e4,
        crack_v = 200,
        A = 0.4 * 0.05
    )
hydrogen_results = hydrogen_simulator.simulate(target_W=0.0, tolerance=1e-3)

# hydrogen2_simulator = IsentropicFlowSimulator(
#         gas_name="hydrogen",
#         P0=8e6,      # 8MPa
#         T0=293,    # 293 K
#         u0=0.0,      # 初始静止
#         dt=1e-4,     # 时间步长
#         max_steps=1e4,
#         crack_v = 200
#     )
# hydrogen2_results = hydrogen2_simulator.simulate(target_W=0.0, tolerance=1e-3)


# methane_simulator = IsentropicFlowSimulator(
#         gas_name="methane",
#         P0=8e6,      # 8MPa
#         T0=293,    # 293 K
#         u0=0.0,      # 初始静止
#         dt=1e-5,     # 时间步长
#         max_steps=1e5
#     )
# methane_results = methane_simulator.simulate(target_W=0.0, tolerance=1e-3)


CrackVelocityCalculator_instance1 = CrackVelocityCalculator(
    t=6.3,        # 管道壁厚 (mm)
    D=463.0,       # 管道直径 (mm)
    YS=430.0,      # 屈服应力 (MPa)
    TS=530.0,      # 抗拉强度 (MPa)
    Pmin=1,         # 最小裂纹尖端扩展压力 (MPa)
    Pmax=8.24,        # 最大裂纹尖端扩展压力 (MPa)
    Cv=23.5,      # 夏比冲击吸收能量 (J)
    Ac=100.0      # 预裂纹夏比试样韧带面积 (mm²)
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

# draw([
#     {"x":hydrogen_a_minus_u, "y":hydrogen_scaled_gas_pressure, "style":'b-', "label":'hydrogen'},
#     {"x":hydrogen2_a_minus_u, "y":hydrogen2_scaled_gas_pressure, "style":'r--', "label":'hydrogen2'}],
#      xlabel='减压速度 u (m/s)',
#      ylabel='压力 P (MPa)',
#      title='压力-速度关系曲线',
#      xlim=(0, 1500),
#      ylim=(0, 12))

# scaled_crack_pressure = [p for p in crack_results['pressure']]
# scaled_crack_pressure2 = [p for p in crack_results2['pressure']]

# draw([
#     {"x":crack_results['velocity'], "y":scaled_crack_pressure, "style":'r--', "label":'crack'},
#     {"x":crack_results2['velocity'], "y":scaled_crack_pressure2, "style":'b-', "label":'crack2'}],
#      xlabel='减压速度 u (m/s)',
#      ylabel='压力 P (MPa)',
#      title='压力-速度关系曲线',
#      xlim=(0, 1500),
#      ylim=(0, 12))

# instance_list1 = [
#     {
#         "id":1,
#         "t": 2.37,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm)
#         "P": 8.25e6    # 管道压力 (Pa)
#     },
#     {
#         "id":2,
#         "t": 2.37,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 8.86e6    # 管道压力 (Pa)
#     },
#     {
#         "id":3,
#         "t": 1.65,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 5.94e6    # 管道压力 (Pa)
#     },
#     {
#         "id":4,
#         "t": 2.37,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 8.84e6    # 管道压力 (Pa)
#     },
#     {
#         "id":5,
#         "t": 1.6,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 9.8e6    # 管道压力 (Pa)
#     },
#     {
#         "id":6,
#         "t": 2.72,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 8.24e6    # 管道压力 (Pa)
#     },
#     {
#         "id":7,
#         "t": 1.6,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 10.08e6    # 管道压力 (Pa)
#     },
#     {
#         "id":8,
#         "t": 1.6,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 10.3e6    # 管道压力 (Pa)
#     }
# ]

# instance_list2 = [
#     {
#         "id":1,
#         "t": 6.3,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm)
#         "P": 8.25e6    # 管道压力 (Pa)
#     },
#     {
#         "id":2,
#         "t": 6.3,        # 管道壁厚 (mm)
#         "D": 462.6,       # 管道直径 (mm),
#         "P": 8.86e6    # 管道压力 (Pa)
#     }
# ]

# for idx, instance in enumerate(instance_list1):
#     hydrogen_simulator = IsentropicFlowSimulator(
#         gas_name="hydrogen",
#         P0=instance["P"],      # 8MPa
#         T0=293,    # 293 K
#         u0=0.0,      # 初始静止
#         dt=1e-4,     # 时间步长
#         max_steps=1e4,
#         crack_v = 200
#     )
#     hydrogen_results = hydrogen_simulator.simulate(target_W=0.0, tolerance=1e-3)
#     crack_calculator = CrackVelocityCalculator(
#         t=instance["t"],
#         D=instance["D"],
#         YS=430,
#         TS=530,
#         Pmin=1,
#         Pmax=12,
#         Cv=200.0,
#         Ac=100.0
#     )
#     crack_results = crack_calculator.calculate()
#     hydrogen_a_minus_u = [a - u for a, u in zip(hydrogen_results['sound_speed'], hydrogen_results['velocity'])]
#     hydrogen_scaled_gas_pressure = [p / 1e6 for p in hydrogen_results['pressure']]
#     # hydrogen2_a_minus_u = [a - u for a, u in zip(hydrogen2_results['sound_speed'], hydrogen2_results['velocity'])]
#     # hydrogen2_scaled_gas_pressure = [p / 1e6 for p in hydrogen2_results['pressure']]
#     # methane_a_minus_u = [a - u for a, u in zip(methane_results['sound_speed'], methane_results['velocity'])]
#     # methane_scaled_gas_pressure = [p / 1e6 for p in methane_results['pressure']]
#     scaled_crack_pressure = [p for p in crack_results['pressure']]
#     draw([{"x":hydrogen_a_minus_u, "y":hydrogen_scaled_gas_pressure, "style":'b-', "label":'hydrogen'},
#         #   {"x":methane_a_minus_u, "y":methane_scaled_gas_pressure, "style":'g-', "label":'methane'},
#         {"x":crack_results['velocity'], "y":scaled_crack_pressure, "style":'r--', "label":'crack'}],
#         xlabel='减压速度 u (m/s)',
#         ylabel='压力 P (MPa)',
#         title=f'压力-速度关系曲线 (实例 {instance["id"]})',
#         xlim=(0, 1500),
#         ylim=(0, 12))
    