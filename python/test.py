from nikfemm import MultiLayerCurrentDensitySimulation

simulation = MultiLayerCurrentDensitySimulation(2, [70e-6, 70e-6])


length = 1
width = 0.001

simulation.draw_rectangle([0, 0], [length, width], 0)
simulation.draw_rectangle([0, 0], [length, width], 1)

simulation.draw_region([0.5 * length, 0.5 * width], 5.95e7, 0)
simulation.draw_region([0.5 * length, 0.5 * width], 5.95e7, 1)

simulation.add_interconnection([length - 0.5 * width, 0.5 * width], [0.5 * width, 0.5 * width], 0, 1, 1)

system = simulation.generate_system(False, 1)

simulation.set_voltage(system, [0, 0.5 * width], -1, 0)
simulation.set_voltage(system, [length, 0.5 * width], 1, 1)

print(system)

simulation.solve(system)

print(simulation.get_layer_voltages(0))