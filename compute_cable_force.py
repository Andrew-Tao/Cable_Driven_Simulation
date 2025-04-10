from cable_functions import get_ALE_elem_in_cable,compute_Q_i_ALE_elem,compute_Q_e_ALE_elem,compute_Q_g_ALE_elem
import numpy as np

q_c = np.array([0,0,0,0,0,0,1,1,0,0,2,2,0,0,3,3,0,0,4,4])
dot_q_c = np.array([0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,2,0])
dot_dot_q_c = np.array([0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,2,0])

def compute_Q_i_individual_cable(q_c,dot_q_c,dot_dot_q_c):
    n_elements = int(len(q_c) / 4)

    Q_i = np.array([])
    for i in range(n_elements - 1):

        if i == 0:
            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_i = compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)[0:4]
            Q_i = np.append(Q_i,individual_Q_i)
        else:
            individual_q_c = get_ALE_elem_in_cable(q_c, i - 1)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i - 1)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i - 1)
            individual_Q_i = compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)
            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_i = individual_Q_i[4 : 8] + compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c,
                                                                         individual_dot_dot_q_c)[0:4]

            Q_i = np.append(Q_i, individual_Q_i)

    individual_q_c = get_ALE_elem_in_cable(q_c, n_elements - 2)
    individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, n_elements - 2)
    individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, n_elements - 2)
    individual_Q_i = compute_Q_i_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)[4:8]
    Q_i = np.append(Q_i, individual_Q_i)
    return Q_i

def compute_Q_e_individual_cable(q_c,dot_q_c,dot_dot_q_c):
    n_elements = int(len(q_c) / 4)

    Q_e = np.array([])
    for i in range(n_elements - 1):

        if i == 0:
            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_e = compute_Q_e_ALE_elem(individual_q_c, individual_dot_q_c, 1,1,1)[0:4]
            Q_e = np.append(Q_e,individual_Q_e)
        else:
            individual_q_c = get_ALE_elem_in_cable(q_c, i - 1)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i - 1)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i - 1)
            individual_Q_e = compute_Q_e_ALE_elem(individual_q_c, individual_dot_q_c, 1,1,1)
            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_e = individual_Q_e[4 : 8] + compute_Q_e_ALE_elem(individual_q_c, individual_dot_q_c, 1,1,1)[0:4]

            Q_e = np.append(Q_e, individual_Q_e)

    individual_q_c = get_ALE_elem_in_cable(q_c, n_elements - 2)
    individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, n_elements - 2)
    individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, n_elements - 2)
    individual_Q_e = compute_Q_e_ALE_elem(individual_q_c, individual_dot_q_c, 1,1,1)[4:8]
    Q_e = np.append(Q_e, individual_Q_e)
    return Q_e

def compute_Q_g_individual_cable(q_c,dot_q_c,dot_dot_q_c):
    n_elements = int(len(q_c) / 4)

    Q_g = np.array([])

    for i in range(n_elements - 1):
        #print("i", i)
        if i == 0:
            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_g = compute_Q_g_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)[0:4]
            Q_g = np.append(Q_g,individual_Q_g)
        else:
            individual_q_c = get_ALE_elem_in_cable(q_c, i - 1)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i - 1)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i - 1)
            individual_Q_g = compute_Q_g_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)
            #print("individual_Q_g", individual_Q_g)
            individual_q_c = get_ALE_elem_in_cable(q_c, i)
            individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, i)
            individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, i)
            individual_Q_g = individual_Q_g[4 : 8] + compute_Q_g_ALE_elem(1, 1, individual_q_c, individual_dot_q_c,
                                                                         individual_dot_dot_q_c)[0:4]
            """
            print("individual_Q_g", compute_Q_g_ALE_elem(1, 1, individual_q_c, individual_dot_q_c,
                                                                         individual_dot_dot_q_c))
            print("individual_Q_g", individual_Q_g)
            """
            Q_g = np.append(Q_g, individual_Q_g)

    individual_q_c = get_ALE_elem_in_cable(q_c, n_elements - 2)
    individual_dot_q_c = get_ALE_elem_in_cable(dot_q_c, n_elements - 2)
    individual_dot_dot_q_c = get_ALE_elem_in_cable(dot_dot_q_c, n_elements - 2)
    individual_Q_g = compute_Q_g_ALE_elem(1, 1, individual_q_c, individual_dot_q_c, individual_dot_dot_q_c)[4:8]
    Q_g = np.append(Q_g, individual_Q_g)
    return Q_g

#print(compute_Q_g_individual_cable(q_c,dot_q_c,dot_dot_q_c))



