import numpy as np
from test import func_i,func_g

p = 4.4
p1 = 2
p2 = 10
dot_p1 = 0
dot_p2 = 0
r1 = np.array([0,2,4])
r2 = np.array([0,0,1])

dot_r1 = np.array([3,4,5])
dot_r2 = np.array([2,4,6])
q_c = np.array([0,0,1,1,0,0,2,2])
dot_q_c = np.array([0,0,1,dot_p1,0,0,2,dot_p2])
dot_dot_q_c = np.array([0,0,1,0,0,0,1,0])
density = 3
cross_section_area = 4
youngs_modulus = 1
damping_coefficient = 0.2
def compute_Q_i_ALE_elem(density,cross_section_area,q_c,dot_q_c,dot_dot_q_c):
    # generate p
    p1 = q_c[3]
    p2 = q_c[7]
    #print("p2", p2, "p1", p1)
    integral_steps = 1000
    integral_length =(p2-p1) / (integral_steps)
    p = np.linspace(p1,p2,integral_steps)
    result = np.zeros(8)
    for i in p:
        result += func_i(i,q_c,dot_q_c,dot_dot_q_c,density,cross_section_area) * integral_length

    return result


def compute_Q_e_ALE_elem(q_c,dot_q_c,cross_section_area,youngs_modulus,damping_coefficient):
    p1 = q_c[3]
    p2 = q_c[7]
    r1 = q_c[0:3]
    r2 = q_c[4:7]

    # compute epsilon
    epsilon = np.linalg.norm(r2-r1) / (p2 - p1) - 1


    # compute k
    if (epsilon > 0): k = 1
    else: k = 0

    # compute B
    distance = np.linalg.norm(r2-r1)
    length = p2-p1

    t1 = -(r2-r1).T/ (distance*length)
    t2 = np.array([ distance / (length**2)])
    t3 = (r2 - r1).T / (distance * length)
    t4 = np.array([ -distance / (length**2)])

    B = np.concatenate((t1,t2,t3,t4))



    # compute dot_epsilon
    dot_epsilon = B @ dot_q_c

    # compute K
    K = k * epsilon - damping_coefficient * dot_epsilon



    result = (p2-p1) * cross_section_area * youngs_modulus * K * B.T # K = (k epsilon - c dot_epsilon)
    return result

def compute_Q_g_ALE_elem(density,cross_section_area,q_c,dot_q_c,dot_dot_q_c):
    # generate p
    p1 = q_c[3]
    p2 = q_c[7]
    integral_steps = 1000
    integral_length =(p2-p1) / (integral_steps)
    p = np.linspace(p1,p2,integral_steps)
    result = np.zeros(8)
    for i in p:
        result += func_g(i,q_c,dot_q_c,dot_dot_q_c,density,cross_section_area) * integral_length

    return result


def get_ALE_elem_in_cable(generalized_coordinate,n):
    return generalized_coordinate[4*n:4*n+8]

