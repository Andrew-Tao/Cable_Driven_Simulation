
from numba.cpython.randomimpl import laplace_impl1
from tqdm.notebook import tqdm_notebook

from cable import ALE_Cable
from cable_functions import compute_Q_i_ALE_elem, dot_dot_q_c, get_ALE_elem_in_cable
from compute_cable_force import compute_Q_i_individual_cable,compute_Q_e_individual_cable,compute_Q_g_individual_cable
from rigid_link import Rigid_Link
import numpy as np

class Simulator():

    def __init__(self,cable_collection,rigid_link_collection,q)->None:
        self.cable_collection = cable_collection
        self.rigid_link_collection = rigid_link_collection
        self.q = q

    @classmethod
    def simulate_system(cls,cable,rigid_link):
        cable_collection = cable
        rigid_link = rigid_link
        cable_1 = cable_collection[0]
        cable_2 = cable_collection[1]
        cable_3 = cable_collection[2]
        q_c_1 = np.array(cable_1.q_c)
        q_c_2 = np.array(cable_2.q_c)
        q_c_3 = np.array(cable_3.q_c)
        r_b = rigid_link.q_b
        q = np.concatenate((q_c_1.ravel(), q_c_2.ravel(), q_c_3.ravel(), r_b.ravel())).T

        #print("q", q)

        return cls(cable,rigid_link,q)

    # End_Point_Constrain
    def constrain_h_1(self):
        constrain_h_1 = np.array([])
        for i in range(42):
            num = 0
            for cable in self.cable_collection.copy():
                if i % 7 == 2 or i % 7 == 3 or i % 7 == 4 or i % 7 == 5:
                    r_i = cable.global_coordinate[:,i]
                    r_bj = self.rigid_link_collection.position.copy()[:,int(i/7)]
                    alpha = (-(2/3)*num - (1/2))* np.pi
                    rho = np.array([(0.05/7)*(i%7)-0.025,0.02*np.sin(alpha),0.02*np.cos(alpha)])
                    A_j =  self.rigid_link_collection.compute_matrix_A(int(i/7))
                    num +=1
                    constrain_h_1_individual = r_i - (r_bj + (A_j @ rho))
                    constrain_h_1 = np.concatenate((constrain_h_1,constrain_h_1_individual))
                    #print("rigid_link",int(i/7),"index",i%7,"num",num)
                    #print("ri",cable.global_coordinate[:,i] -(r_bj+(A_j @ rho)) )
                    #print("r_bj",r_bj+(A_j @ rho))
                else:
                    constrain_h_1 = np.concatenate((constrain_h_1, np.array([0.0,0.0,0.0])))
    def constrain_h_2(self):
        constrain_h_2_tem = np.array([])

        for i in range(42):
            num = 0
            for cable in self.cable_collection.copy():
                if i % 7 == 2 or i % 7 == 3 or i % 7 == 4 or i % 7 == 5:
                    r_i = cable.global_coordinate[:, i]
                    r_bj = self.rigid_link_collection.position.copy()[:, int(i / 7)]
                    alpha = (-(2 / 3) * num - (1 / 2)) * np.pi
                    rho = np.array([(0.05 / 7) * (i % 7) - 0.025, 0.02 * np.sin(alpha),0.02 * np.cos(alpha) ])
                    A_j = self.rigid_link_collection.compute_matrix_A(int(i / 7))
                    rho_n = np.array([1,0,0])
                    num += 1
                    constrain_h_2_individual_tem = np.dot((r_i - (r_bj + (A_j @ rho))),(A_j @ rho_n))
                    constrain_h_2_individual = np.full((3), constrain_h_2_individual_tem)
                    constrain_h_2_tem = np.concatenate((constrain_h_2_tem,constrain_h_2_individual))
                    # print("rigid_link",int(i/7),"index",i%7,"num",num)
                    # print("ri",cable.global_coordinate[:,i] )
                    # print("r_bj",r_bj+(A_j @ rho))
                else:
                    constrain_h_2_tem = np.concatenate((constrain_h_2_tem,np.array([0.0,0.0,0.0])))

        constrain_h_2 = np.zeros_like(constrain_h_2_tem)
        for i in range(42):
            constrain_h_2[i] = constrain_h_2_tem[i*3]
            constrain_h_2[i+42] = constrain_h_2_tem[i*3+1]
            constrain_h_2[i+84] = constrain_h_2_tem[i*3+2]
        cable1 = self.cable_collection[0]
        print("debug3",cable1.material_coordinate)

        return(constrain_h_2)

    #Pull Cable_1, Free Another Two cable
    def constrain_of_material_flow(self):
        #Fix Tip Material Flow
        constrain_p_tem = np.array([])
        """
        for cable in self.cable_collection.copy():
            cable.material_cooridnate[-1] = 0.05
            cable.generalized_coordinate[-1] = 0.05
        """
        target_material_cooridnate = np.array([])
        #cable1 = self.cable_collection[2]
        #print(cable1.material_coordinate)
        for cable in self.cable_collection.copy():
            print("debug0",cable.material_coordinate)
            target_material_cooridnate = np.concatenate((target_material_cooridnate,cable.material_coordinate.copy()))
            constrain_p_tem = np.concatenate((constrain_p_tem,cable.material_coordinate.copy()))


        constrain_p = target_material_cooridnate - constrain_p_tem


        return(constrain_p)
    def constrain_A (self):
        constrain_A = np.array([])
        link = self.rigid_link_collection
        for i in range(6):
            constrain_A_ind = np.dot(link.orientation.copy()[:,i],link.orientation.copy()[:,i]) - 1
            constrain_A_tem = np.full((4), constrain_A_ind)
            constrain_A = np.concatenate((constrain_A,constrain_A_tem))

        print(constrain_A.shape)
        return(constrain_A)

    # I am not sure whether this constrain is correct
    def constrain_J(self):

        constrain_J = np.array([])

        rigid_link = self.rigid_link_collection
       # print(rigid_link.compute_matrix_A(0))
       # print(rigid_link.compute_matrix_A(1))
       # print(rigid_link.compute_matrix_A(2))
        A_j = rigid_link.compute_matrix_A(0)
        end_point_position = np.array([0,0,0]) # Need Fix, you are supposed to get the data from cable
      #  print(end_point_position)
        def norm(vector):
            vector = vector/np.linalg.norm(vector)
            return vector

        constrain_J_ind = end_point_position + (0.025*(norm(A_j @ np.array([1.0,0.0,0.0]))))
       # print("ind",A_j @ np.array([1.0,0.0,0.0]))
       # print(rigid_link.position[:,0].copy())
        constrain_J = np.concatenate((constrain_J,constrain_J_ind-rigid_link.position[:,0].copy()))
        for i in range(1,6):
            #print(i)
            A_j_last = rigid_link.compute_matrix_A(i-1)
            A_j = rigid_link.compute_matrix_A(i)
            constrain_J_ind = constrain_J_ind + (0.025*(norm(A_j_last @ np.array([1.0,0.0,0.0])))) + (0.025*(norm(A_j @ np.array([1.0,0.0,0.0]))))
            constrain_J = np.concatenate((constrain_J, constrain_J_ind-rigid_link.position[:,i].copy()))
            #print("ind", A_j @ np.array([1.0,0.0,0.0]))
           # print(rigid_link.position[:, i].copy())


        #print("J",constrain_J, constrain_J.shape)
        return(constrain_J)


    def contact_between_cable_and_hole(self):
            Q_nc_tem= np.array([])
            Q_nh= np.array([])
            for i in range(42):
                num = 0
                for cable in self.cable_collection.copy():
                    if i % 7 == 2 or i % 7 == 3 or i % 7 == 4 or i % 7 == 5:
                        r_c = cable.global_coordinate[:, i]
                        #print("debug2",cable.global_coordinate)
                        r_bj = self.rigid_link_collection.position.copy()[:, int(i / 7)]
                        orientation = self.rigid_link_collection.orientation.copy()
                        lam0 = orientation[0, int(i / 7)]
                        lam1 = orientation[1, int(i / 7)]
                        lam2 = orientation[2, int(i / 7)]
                        lam3 = orientation[3, int(i / 7)]

                        R_h = 0.002
                        R_c = 0.0015
                        alpha = (-(2 / 3) * num - (1 / 2)) * np.pi
                        rho = np.array([(0.05 / 7) * (i % 7) - 0.025, 0.02 * np.sin(alpha), 0.02 * np.cos(alpha)])

                        A_j = self.rigid_link_collection.compute_matrix_A(int(i / 7))

                        r_h = r_bj + (A_j @ rho)

                        #print(i,r_h - r_c)
                        delta = np.linalg.norm(r_h-r_c) - (R_h - R_c)
                        K = 1e6
                        C = 10e6
                        dot_delta = 0
                        if delta < 0:      # This part need fix
                            F_n= 0
                        else:
                            F_n = ((1/2) * K * (R_c /100)) + (K * delta) + C * dot_delta

                        n_c = (r_c - r_h) / np.linalg.norm(r_h - r_c)
                        n_h = (r_h - r_c) / np.linalg.norm(r_h - r_c)

                        A_partial_lam0 = np.array([
                            [4*lam0 , -2*lam3 , 2*lam2],
                            [2*lam3, 4*lam0, -2*lam1],
                            [-2* lam2, 2* lam1, 4*lam0]
                        ])
                        A_partial_lam1 = np.array([
                            [4*lam1, 2*lam2, 2*lam3],
                            [2*lam2, 0, -2*lam0],
                            [2*lam3, 2*lam0,0]
                        ])
                        A_partial_lam2 = np.array([
                            [0,2*lam1,2*lam0],
                            [2*lam1,4*lam2,2*lam3],
                            [-2*lam0,2*lam3,0]
                        ])
                        A_partial_lam3 = np.array([
                            [0,-2*lam0,2*lam1],
                            [2*lam0,0,2*lam2],
                            [2*lam1,2*lam2,4*lam3]
                        ])

                        N_c = np.concatenate((np.identity(3), np.zeros((3,1))), axis=1)
                        #print(A_partial_lam0 @ rho)
                        t_0_tem = (A_partial_lam0 @ rho)
                        t_0 = np.array([[t_0_tem[0]],[t_0_tem[1]],[t_0_tem[2]]])
                        t_1_tem = (A_partial_lam1 @ rho)
                        t_1 = np.array([[t_1_tem[0]],[t_1_tem[1]],[t_1_tem[2]]])
                        t_2_tem = (A_partial_lam2 @ rho)
                        t_2 = np.array([[t_2_tem[0]],[t_2_tem[1]],[t_2_tem[2]]])
                        t_3_tem = (A_partial_lam3 @ rho)
                        t_3 = np.array([[t_3_tem[0]], [t_3_tem[1]], [t_3_tem[2]]])


                        N_h = np.concatenate(
                        (np.identity(3),
                                t_0,
                                t_1,
                                t_2 ,
                                t_3 ),axis=1)
                        Q_nc_individual = F_n * N_c.T @ n_c
                        Q_nh_individual = F_n * N_h.T @ n_h
                        #print(Q_nh_individual)

                        Q_nc_tem = np.concatenate((Q_nc_tem, Q_nc_individual))
                        Q_nh = np.concatenate((Q_nh, Q_nh_individual))


                        num += 1
                    else:
                        Q_nc_tem = np.concatenate((Q_nc_tem, np.zeros((4))))
                        Q_nh = np.concatenate((Q_nh, np.zeros((7))))


                # Resort Q_nc based on cable_order
                Q_nc = np.zeros_like(Q_nc_tem)
                #print(Q_nc)
            #print("debug",Q_nc_tem.shape)

                
            for i in range(42):
                Q_nc[i] = Q_nc_tem[i*3]
                Q_nc[i+42] = Q_nc_tem[i*3+1]
                Q_nc[i+84] = Q_nc_tem[i*3+2]

            def compute_Q_nb_from_Q_nh(Q_nh):

                Q_nb_link_1 = np.zeros((7))
                Q_nb_link_2 = np.zeros((7))
                Q_nb_link_3 = np.zeros((7))
                Q_nb_link_4 = np.zeros((7))
                Q_nb_link_5 = np.zeros((7))
                Q_nb_link_6 = np.zeros((7))

                for i in range(21):
                    Q_nb_link_1 += Q_nh[i * 7:i * 7 + 7]
                    Q_nb_link_2 += Q_nh[i * 7 + 147:i * 7 + 7 + 147]
                    Q_nb_link_3 += Q_nh[i * 7 + (147 * 2):i * 7 + 7 + (147 * 2)]
                    Q_nb_link_4 += Q_nh[i * 7 + (147 * 3):i * 7 + 7 + (147 * 3)]
                    Q_nb_link_5 += Q_nh[i * 7 + (147 * 4):i * 7 + 7 + (147 * 4)]
                    Q_nb_link_6 += Q_nh[i * 7 + (147 * 5):i * 7 + 7 + (147 * 5)]
                Q_nb = np.concatenate((
                    Q_nb_link_1,
                    Q_nb_link_2,
                    Q_nb_link_3,
                    Q_nb_link_4,
                    Q_nb_link_5,
                    Q_nb_link_6))

                return Q_nb
            Q_nb = compute_Q_nb_from_Q_nh(Q_nh)

            return {"Q_nc":Q_nc,"Q_nh":Q_nh,"Q_nb":Q_nb}




    def run_simulation(cls,dt,time):
        n = 0,
        t = 0,

        q = cls.q
        dot_q = np.zeros_like(q) #Need More Fix, The initial condition
        dot_dot_q = np.zeros_like(q)  # Need more fix, the initial condition
        q_history_collection = np.array([q])
        dot_q_history_collection = np.array([dot_q])
        dot_dot_q_history_collection = np.array([dot_dot_q])
        Lagrange_coeff = 1.0

        # compute dot_q_c & dot_dot_q_c from data collection/ estimate q_c dot q_c based on qn+1

        dot_q = (q - q_history_collection[-1]) / dt
        dot_dot_q = (dot_q - (dot_q_history_collection[-1])) / dt

        # compute residual vector R
        #print(len(q))
        q_c_1 = q[0:168]
        dot_q_c_1 = dot_q[0:168]

        dot_dot_q_c_1 = dot_dot_q[0:168]

        q_c_2 = q[168:336]
        dot_q_c_2 = dot_q[168:336]
        dot_dot_q_c_2 = dot_dot_q[168:336]
        q_c_3 = q[336:504]
        dot_q_c_3 = dot_q[336:504]
        dot_dot_q_c_3 = dot_dot_q[336:504]
        """

        for i in range(len(q_c_2)):
            if (i%4)==3:
                print(q_c_2[i])
        """
        Q_i_1 = compute_Q_i_individual_cable(q_c_1,dot_q_c_1,dot_dot_q_c_1)
        Q_i_2 = compute_Q_i_individual_cable(q_c_2,dot_q_c_2,dot_dot_q_c_2)
        Q_i_3 = compute_Q_i_individual_cable(q_c_3,dot_q_c_3,dot_dot_q_c_3)
        Q_i = np.concatenate((Q_i_1, Q_i_2, Q_i_3))
        Q_e_1 = compute_Q_e_individual_cable(q_c_1,dot_q_c_1,dot_dot_q_c_1)
        Q_e_2 = compute_Q_e_individual_cable(q_c_2,dot_q_c_2,dot_dot_q_c_2)
        Q_e_3 = compute_Q_e_individual_cable(q_c_3,dot_q_c_3,dot_dot_q_c_3)
        Q_e = np.concatenate((Q_e_1, Q_e_2, Q_e_3))
        Q_g_1 = compute_Q_g_individual_cable(q_c_1,dot_q_c_1,dot_dot_q_c_1)
        Q_g_2 = compute_Q_g_individual_cable(q_c_2,dot_q_c_2,dot_dot_q_c_2)
        Q_g_3 = compute_Q_g_individual_cable(q_c_3,dot_q_c_3,dot_dot_q_c_3)
        Q_g = np.concatenate((Q_g_1, Q_g_2, Q_g_3))
        Contact_force = cls.contact_between_cable_and_hole()
        Q_nc = Contact_force["Q_nc"]
        Q_nh = Contact_force["Q_nh"]
        Q_nb = Contact_force["Q_nb"]

        rigid_link = cls.rigid_link_collection
        M = rigid_link.compute_M_generalized()
        dot_dot_q_b = rigid_link.dot_dot_q_b
        Q_I = rigid_link.compute_Q_I_generalized()
        Q_b = rigid_link.compute_Q_B_generalized()
        H = rigid_link.compute_H_generalized()

        Q_dynamic = np.concatenate((Q_i + Q_e, M @ dot_dot_q_b + H.T @ Q_I))
        Q_ext = np.concatenate((Q_g + Q_nc, H.T @ Q_I + Q_nb))
        constrain_h = cls.constrain_h_2()
        constrain_p = cls.constrain_of_material_flow()
        constrain_A = cls.constrain_A()
        constrain_J = cls.constrain_J()



        # Compute Constrain_Jacob ???????

        constrain_jacob =3.0
        constrain = np.concatenate((constrain_J, constrain_A, constrain_h, constrain_p))
        F = Q_dynamic + Q_ext +  Lagrange_coeff * constrain_jacob
        Q = constrain

        R = np.array([F,Q])

        print("constrain",constrain.shape)


        #print(Q_constrain.shape)
        """
        print(H.shape)
        print(M.shape)
        print(dot_dot_q_b.shape)
        print(Q_I.shape)
        print(Q_nh.shape)

        print(Q_nb.shape)
        print(Q_dynamic.shape)
        print(Q_ext.shape)
        """
        print(Q_dynamic.shape)
        print(H.shape)



        return True







"""
def simulate (cable,rigid_link,dt,time):



    q_combine = np.concatenate((
        cable_collection[0].generalized_coordinate,
        cable_collection[1].generalized_coordinate,
        cable_collection[2].generalized_coordinate,
        rigid_link_1.q_b,
    )).T

    data_collection = {"time":np.array([0.0]),"q_combine_collection":np.array([q_combine])}

    print(q_combine)
   
    
    total_step = int (time / dt)
    iteration = 0
    
    
    
    for t in range(total_step):



def compute_Q_i_individual_cable(cable, data_collection, dt):
    n_elements = cable.n_elem.copy()
    generalized_coordinate = cable.generalized_coordinate.copy()
    Q_i = np.array([])

    # compute dot_q_c & dot_dot_q_c from data collection/ estimate q_c dot q_c based on qn+1
    q_c = generalized_coordinate
    dot_q_c = (generalized_coordinate - data_collection[-1]) / dt
    dot_dot_q_c = (dot_q_c - ((data_collection[-1] - data_collection[-2]) / dt)) / dt

    for i in range(n_elements - 1):
        # wrong: write it in integral form
        if i == 0:
            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_i = compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)

            Q_i.append(individual_Q_i)
        else:
            individual_q_c = get_ALE_elem_in_cable(q_c, i - 1)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i - 1)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i - 1)
            individual_Q_i = compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)

            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_i = individual_Q_i[4, 8] + compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c,
                                                                         individual_dot_dot_q_c)[0:4]

            Q_i.append(individual_Q_i)

        individual_q_c = get_ALE_elem_in_cable(q_c, n_elements - 1)
        individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, n_elements - 1)
        individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, n_elements - 1)
        individual_Q_i = compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)[4:8]
        Q_i.append(individual_Q_i)
    return Q_i


#simulate(cable_collection,rigid_link_1,dt=1,time = 100)

# data_collection (time_step,(cable)) 2-Dimensional
"""

