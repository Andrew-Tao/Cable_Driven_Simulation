from cable import ALE_Cable
from Simulator import Simulator
import numpy as np

from cable_functions import youngs_modulus
from rigid_link import Rigid_Link
youngs_modulus = 1.39 * (10 **9)
cross_section_area_cable = np.pi * (10 ** (-8))
direction = np.array([0,0,-1])
length = 0.3

cable_1 = ALE_Cable.normal_cable(
    42,
    cross_section_area_cable,
    youngs_modulus,
    end_point_position = np.array([0,0.02,0]),
    length = length,
    director = direction,
)
cable_2 = ALE_Cable.normal_cable(
    42,
    cross_section_area_cable,
    youngs_modulus,
    end_point_position = np.array([0.0173205,-0.01,0]),
    length = length,
    director = direction,)
cable_3 = ALE_Cable.normal_cable(
    42,
    cross_section_area_cable,
    youngs_modulus,
    end_point_position = np.array([-0.0173205,-0.01,0]),
    length = length,
    director = direction,)

cable_collection = [cable_1, cable_2, cable_3]

#print("debug4",cable_1.material_coordinate)

rigid_link_collection = Rigid_Link.rigid_link(6,length,direction,np.array([0,0,-0.025]))

cable_driven_snake = Simulator.simulate_system(cable_collection,rigid_link_collection)
#print(cable_driven_snake)

print(cable_driven_snake.run_simulation(0.01,3))

