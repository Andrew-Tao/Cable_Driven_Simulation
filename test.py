import numpy as np

p = 4.4
p1 = 2
p2 = 10
dot_p1 = 1
dot_p2 = 1
r1 = np.array([0,2,4])
r2 = np.array([0,0,1])

dot_r1 = np.array([3,4,5])
dot_r2 = np.array([2,4,6])
q_c = np.array([0,2,4,2,0,0,1,10])
dot_q_c = np.array([3,4,5,dot_p1,2,4,6,dot_p2])


def func_i(
        p: float,
        q_c: np.ndarray,
        dot_q_c: np.ndarray,
        dot_dot_q_c: np.ndarray,
        rho,
        A,
):

    assert q_c.shape == (8,)
    assert dot_q_c.shape == (8,)
    assert dot_dot_q_c.shape == (8,)

    p1 = q_c[3]
    p2 = q_c[7]
    r1 = q_c[0:3]
    r2 = q_c[4:7]
    dot_p1 = dot_q_c[3]
    dot_p2 = dot_q_c[7]
    dot_r1 = dot_q_c[0:3]
    dot_r2 = dot_q_c[4:7]

    # Compute Xi
    Xi = (p - p1) / (p2 - p1)
    N1 = 1 - Xi
    N2 = Xi
    r = N1 * r1 + N2 * r2
    N1_partial_p1 = (-p + p2) / ((p2 - p1) ** 2)
    N1_partial_p2 = (p - p1) / ((p2 - p1) ** 2)
    N2_partial_p1 = (p2 - p) / ((p2 - p1) ** 2)
    N2_partial_p2 = (p1 - p) / ((p2 - p1) ** 2)

    t1 = N1 * np.identity(3)
    t2 = np.array([ N1_partial_p1 * r1 + N2_partial_p1 * r2 ])
    t3 = N2 * np.identity(3)
    t4 = np.array([N1_partial_p2 * r1 + N2_partial_p2 * r2])

    N = np.concatenate((t1,t2.T,t3,t4.T), axis=1)
    a = (N1_partial_p2 * dot_p1 + N1_partial_p2 * dot_p2) * dot_r1 + (
                N2_partial_p2 * dot_p1 + N2_partial_p2 * dot_p2) * dot_r2

    #dot_r = N @ dot_q_c
    dot_dot_r = N @ dot_dot_q_c + a
    result = (rho * A) * N.T @ dot_dot_r

    return result

def func_g(
        p: float,
        q_c: np.ndarray,
        dot_q_c: np.ndarray,
        dot_dot_q_c: np.ndarray,
        rho,
        A,
        gravity_acceleration = -9.8
):

    assert q_c.shape == (8,)
    assert dot_q_c.shape == (8,)
    assert dot_dot_q_c.shape == (8,)

    p1 = q_c[3]
    p2 = q_c[7]
    r1 = q_c[0:3]
    r2 = q_c[4:7]
    dot_p1 = dot_q_c[3]
    dot_p2 = dot_q_c[7]
    dot_r1 = dot_q_c[0:3]
    dot_r2 = dot_q_c[4:7]

    # Compute Xi
    Xi = (p - p1) / (p2 - p1)
    N1 = 1 - Xi
    N2 = Xi
    r = N1 * r1 + N2 * r2
    N1_partial_p1 = (-p + p2) / ((p2 - p1) ** 2)
    N1_partial_p2 = (p - p1) / ((p2 - p1) ** 2)
    N2_partial_p1 = (p2 - p) / ((p2 - p1) ** 2)
    N2_partial_p2 = (p1 - p) / ((p2 - p1) ** 2)

    t1 = N1 * np.identity(3)
    t2 = np.array([ N1_partial_p1 * r1 + N2_partial_p1 * r2 ])
    t3 = N2 * np.identity(3)
    t4 = np.array([N1_partial_p2 * r1 + N2_partial_p2 * r2])

    N = np.concatenate((t1,t2.T,t3,t4.T), axis=1)
    a = (N1_partial_p2 * dot_p1 + N1_partial_p2 * dot_p2) * dot_r1 + (
                N2_partial_p2 * dot_p1 + N2_partial_p2 * dot_p2) * dot_r2

    g = np.array([0,0,gravity_acceleration])
    result = (rho * A) * N.T @ g

    return result

#print(func_i(p,q_c,dot_q_c,np.array([1,2,3,4,5,6,7,8]),3,4))