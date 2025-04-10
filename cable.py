import numpy as np
from fontTools.misc.bezierTools import epsilon
from numba.cuda import external_stream
from numpy.typing import NDArray
#class ale_element:

#def compute_internal_force(ALE_Cable):


def allocate_cable (
        n_elements: int,
        cross_section_area : float,
        young_modulus: float,
        end_point_position = np.array([0,0,0]),
        length = 10,
        director = np.array([1,0,0]),
) -> tuple[
    int,
    float,
    float,
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    director = director / np.linalg.norm(director)
    length_elem = length/n_elements
    n_elem = n_elements
    global_coordinate = np.zeros((3,n_elem))
    material_coordinate = np.zeros((n_elem))
    epsilon = np.zeros ((n_elem))
    N_1 = np.zeros((n_elem))
    N_2 = np.zeros((n_elem))
    epsilon[0] = 0
    for i in range (n_elem):
        global_coordinate[:,i] = i*length_elem*director + end_point_position
        material_coordinate[i] = i * length_elem

    generalized_coordinate = np.zeros((n_elem*4))

    for i in range (n_elem*4):
        if i % 4 == 3:
            generalized_coordinate[i] = material_coordinate[int((i-3)/4)]

        else:
            generalized_coordinate[i] = global_coordinate[i%4,np.floor(i/4).astype(int)]


    generalized_coordinate = generalized_coordinate.T
    for i in range (n_elem-1):
        epsilon[i+1] = np.linalg.norm(global_coordinate[:,i+1] - global_coordinate[:,i]) / (material_coordinate[i+1] - material_coordinate[i])

    dot_global_coordinate = np.zeros((3,n_elem))
    dot_dot_global_coordinate = np.zeros((3,n_elem))
    inertia_force = np.zeros((3,n_elem))
    external_force = np.zeros((3,n_elem))
    gravity_force = np.zeros((3,n_elem))
    q_c = generalized_coordinate
    dot_q_c = np.zeros((3,n_elem)) # Need more fix
    dot_dot_q_c = np.zeros((3,n_elem)) # Need more fix
    #print("debug5",material_coordinate)

    return(
        n_elem,
        cross_section_area,
        young_modulus,
        generalized_coordinate,
        global_coordinate,
        dot_global_coordinate,
        dot_dot_global_coordinate,
        material_coordinate,
        epsilon,
        inertia_force,
        external_force,
        gravity_force,
        q_c,
        dot_q_c,
        dot_dot_q_c,
    )
"""
allocate_cable(
    n_elements=18,
    cross_section_area=0.1,
    young_modulus=0.1,
)
"""



class ALE_Cable:
    def __init__ (
            self,
            n_elem: int,
            cross_section_area: float,
            young_modulus: float,
            generalized_coordinate: NDArray[np.float64],
            global_coordinate: NDArray[np.float64],
            dot_global_coordinate: NDArray[np.float64],
            dot_dot_global_coordinate: NDArray[np.float64],
            material_coordinate: NDArray[np.float64],
            epsilon: NDArray[np.float64],    #strain
            inertia_force: NDArray[np.float64],
            external_force: NDArray[np.float64],
            gravity_force: NDArray[np.float64],
            q_c: NDArray[np.float64],
            dot_q_c: NDArray[np.float64],
            dot_dot_q_c: NDArray[np.float64],

    )->None:

        self.n_elem = n_elem
        self.cross_section_area = cross_section_area
        self.young_modulus = young_modulus
        self.generalized_coordinate = generalized_coordinate
        self.global_coordinate = global_coordinate
        self.material_coordinate = material_coordinate
        self.epsilon = epsilon
        self.inertia_force = inertia_force
        self.exteral_force = external_force
        self.gravity_force = gravity_force
        self.dot_global_coordinate = dot_global_coordinate,
        self.dot_dot_global_coordinate = dot_dot_global_coordinate,
        self.q_c = q_c,
        self.dot_q_c = dot_q_c,
        self.dot_dot_q_c = dot_dot_q_c
    @classmethod
    def normal_cable(
            cls,
            n_elements: int,
            cross_section_area: float,
            young_modulus:float,
            end_point_position: NDArray[np.float64],
            length:float,
            director:NDArray[np.float64],
    ) :
        (
            n_elem,
            cross_section_area,
            young_modulus,
            generalized_coordinate,
            global_coordinate,
            dot_global_coordinate,
            dot_dot_global_coordinate,
            material_coordinate,
            epsilon,
            inertia_force,
            external_force,
            gravity_force,
            q_c,
            dot_q_c,
            dot_dot_q_c,
        ) = allocate_cable(
            n_elements,
            cross_section_area,
            young_modulus,
            end_point_position,
            length,
            director,
        )
        return cls(
            n_elem,
            cross_section_area,
            young_modulus,
            generalized_coordinate,
            global_coordinate,
            dot_global_coordinate,
            dot_dot_global_coordinate,
            material_coordinate,
            epsilon,
            inertia_force,
            external_force,
            gravity_force,
            q_c,
            dot_q_c,
            dot_dot_q_c,
        )


#A= ALE_Cable.normal_cable(n_elements=18,cross_section_area=3,young_modulus=0.4)
#print(A.epsilon)