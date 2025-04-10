import numpy as np
from fontTools.misc.bezierTools import epsilon
from mpl_toolkits.mplot3d.proj3d import rotation_about_vector
from numba.cuda import external_stream
from numpy.typing import NDArray
#class ale_element:

#def compute_internal_force(ALE_Cable):

def allocate_rigid_link (
        n_elements: int,
        length: float,
        director: NDArray[np.float64],
        end_point_position: NDArray[np.float64],
) -> tuple[
    int,
    NDArray[np.float64],
    NDArray[np.float64],
]:
    n_elem = n_elements
    position = np.full ((3,n_elem), 3,dtype=np.float64)
    orientation = np.full ((4,n_elem),2, dtype=np.float64)

    for i in range(n_elem):
        orientation[:,i] = np.array([0.0, -0.7071067811865475, 0.0, 0.7071067811865476])

    length_elem = length / n_elem

    for i in range (n_elem):
        position[:,i] = i* length_elem * director + end_point_position
        #print(position)


    q_b = np.array([])
    for i in range (n_elem):
        q_b = np.append(q_b, position[:,i ].T)  # Append position
        q_b = np.append(q_b, orientation[:, i].T)  # Append orientation
    q_b = q_b.T

    dot_q_b = np.zeros_like(q_b)
    dot_q_b[0:4] = 3
    dot_dot_q_b = np.zeros_like(dot_q_b)
    rotational_inertia_matrix = np.identity(3)
    np.fill_diagonal(rotational_inertia_matrix, [1, 1, 2])

    return(
        n_elem,
        position,
        orientation,
        q_b,
        dot_q_b,
        dot_dot_q_b,
        rotational_inertia_matrix,
    )



class Rigid_Link:
    def __init__ (
            self,
            n_elem: int,
            position: NDArray[np.float64],
            orientation: NDArray[np.float64],
            q_b: NDArray[np.float64],
            dot_q_b: NDArray[np.float64],
            dot_dot_q_b: NDArray[np.float64],
            rotational_inertia_matrix: NDArray[np.float64],
    )->None:

        self.n_elem = n_elem
        self.position = position
        self.orientation = orientation
        self.q_b = q_b
        self.dot_q_b = dot_q_b
        self.dot_dot_q_b = dot_dot_q_b
        self.rotational_inertia_matrix = rotational_inertia_matrix

    @classmethod
    def rigid_link(
            cls,
            n_elements: int,
            length: float,
            director: NDArray[np.float64],
            end_point_position: NDArray[np.float64],
    ) :
        (
        n_elem,
        position,
        orientation,
        q_b,
        dot_q_b,
        dot_dot_q_b,
        rotational_inertia_matrix

        ) = allocate_rigid_link(
            n_elements,
            length,
            director,
            end_point_position,
        )
        return cls(
            n_elem,
            position,
            orientation,
            q_b,
            dot_q_b,
            dot_dot_q_b,
            rotational_inertia_matrix,
        )

    def compute_L(self,n_elem: int) -> NDArray[np.float64]:
        orientation = self.orientation.copy()
        lam0 = orientation[0,n_elem]
        lam1 = orientation[1,n_elem]
        lam2 = orientation[2,n_elem]
        lam3 = orientation[3,n_elem]
        L = np.array([
            [-lam1,lam0,lam3,-lam2],
            [-lam2,-lam3,lam0,lam1],
            [-lam3,lam2,-lam1,lam0],
        ])
        return L

    def compute_dot_L(self,n_elem: int) -> NDArray[np.float64]:
        orientation = self.dot_q_b.copy()[7*n_elem+3:7*n_elem+7]
        lam0 = orientation[0]
        lam1 = orientation[1]
        lam2 = orientation[2]
        lam3 = orientation[3]
        dot_L = np.array([
            [-lam1, lam0, lam3, -lam2],
            [-lam2, -lam3, lam0, lam1],
            [-lam3, lam2, -lam1, lam0],
        ])
        return dot_L
    def compute_matrix_A(self,n_elem: int) -> NDArray[np.float64]:
        orientation = self.q_b.copy()[7 * n_elem + 3:7 * n_elem + 7]
        #print(orientation)
        lam0 = orientation[0]
        lam1 = orientation[1]
        lam2 = orientation[2]
        lam3 = orientation[3]
        matrix_A = np.array([
            [2*(lam0**2 + lam1**2)-1 , 2*(lam1*lam2 - lam0*lam3), 2*(lam1*lam3+lam0*lam2)],
            [2*(lam1*lam2+lam0*lam3), 2*(lam0**2+lam2**2)-1,2*(lam2*lam3-lam0*lam1)],
            [2*(lam1*lam3-lam0*lam2), 2*(lam2*lam3 + lam0*lam1),2*(lam0**2+lam3**2)-1]

        ])
        return matrix_A

    def get_dot_orientation(self,n_elem: int) -> NDArray[np.float64]:
        return self.dot_q_b.copy()[7*n_elem+3:7*n_elem+7]

    def compute_H(self,n_elem: int) -> NDArray[np.float64]:
        L = self.compute_L(n_elem)
        H = np.zeros((6,7))
        H[0:3,0:3] = np.identity(3)
        H[3:6,3:7] = L
        return H

    def compute_Q_I (self,n_elem: int) -> NDArray[np.float64]:

        L =self.compute_L(n_elem)
        rotational_inertia_matrix = self.rotational_inertia_matrix.copy()
        #orientation = self.orientation.copy()[:,n_elem]
        dot_L = self.compute_dot_L(n_elem)
        dot_orientation = self.get_dot_orientation(n_elem)
        #print("dot_orient",dot_orientation)
        Q_I = np.zeros((6))
        #print("L",L)
        #print("dot_L",dot_L)
        #print("L @ dot_orientation",dot_L.T @rotational_inertia_matrix @  L @ dot_orientation)
        Q_I[3:6] = 4 * L @ dot_L.T @ rotational_inertia_matrix @ L @ dot_orientation
        return Q_I

    def compute_Z (self,n_elem: int) -> NDArray[np.float64]:
        dot_L = self.compute_dot_L(n_elem)
        Z = np.zeros((6, 6))
        Z[0:3, 0:3] = np.identity(3)
        Z[3:6, 3:6] = self.rotational_inertia_matrix.copy() # Rotatinoal Matrix Need Fix
        return Z

    def compute_M(self,n_elem: int) -> NDArray[np.float64]:
        M = self.compute_H(n_elem).T @ self.compute_Z(n_elem) @ self.compute_H(n_elem)
        return M
    def get_Q_B(self,n_elem: int) -> NDArray[np.float64]:
        Q_b = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        return Q_b
    def compute_M_generalized(self):
        M = np.zeros((42,42))
        M1 = self.compute_M(0)
        M2 = self.compute_M(1)
        M3 = self.compute_M(2)
        M4 = self.compute_M(3)
        M5 = self.compute_M(4)
        M6 = self.compute_M(5)
        M_collection = [M1,M2,M3,M4,M5,M6]
        for i in range(6):
            M[7*i:7*i+7,7*i:7*i+7] = M_collection[i]
        return M
    def compute_H_generalized(self):
        H = np.zeros((36,42))

        H1 = self.compute_H(0)
        H2 = self.compute_H(1)
        H3 = self.compute_H(2)
        H4 = self.compute_H(3)
        H5 = self.compute_H(4)
        H6 = self.compute_H(5)
        H_collection = [H1,H2,H3,H4,H5,H6]
        for i in range(6):
            H[6*i:6*i+6,7*i:7*i+7] = H_collection[i]
        return H
    def compute_Q_I_generalized (self):
        Q1 = self.compute_Q_I(0)
        Q2 = self.compute_Q_I(1)
        Q3 = self.compute_Q_I(2)
        Q4 = self.compute_Q_I(3)
        Q5 = self.compute_Q_I(4)
        Q6 = self.compute_Q_I(5)
        Q_I = np.concatenate((Q1,Q2,Q3,Q4,Q5,Q6))

        return Q_I
    def compute_Q_B_generalized (self):
        Q1 = self.get_Q_B(0)
        Q2 = self.get_Q_B(1)
        Q3 = self.get_Q_B(2)
        Q4 = self.get_Q_B(3)
        Q5 = self.get_Q_B(4)
        Q6 = self.get_Q_B(5)
        Q_B = np.concatenate((Q1,Q2,Q3,Q4,Q5,Q6))

        return Q_B






