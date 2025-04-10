import numpy as np
def length_ray(a,b,k):
    length_ray = a*np.exp((k*b*np.pi)/6)*(np.exp(2*b*np.pi)-1)
    return length_ray

def length_paral(a,b,k):
    length_paral = a*np.exp((k*b*np.pi)/6)*np.sqrt(np.exp((b*np.pi)/3)+1-(np.sqrt(3)*np.exp((b*np.pi)/6)))
    return length_paral

def block(a,b,k): # block provides the four side_length of each half-unit
    individual_unit=np.zeros(4)
    individual_unit[0]=length_ray(a,b,k)
    individual_unit[1]=length_ray(a,b,k+1)
    individual_unit[2]=length_paral(a,b,k)
    individual_unit[3]=length_paral(a,b,k+12)
    return individual_unit
def polar2xy(position):
    x = position[1]* np.cos(position[0])
    y = position[1]* np.sin(position[0])
    position_xy = np.array([x,y])
    return position_xy
def polar_distance(position_A, position_B):
    positionA_xy = polar2xy(position_A)
    positionB_xy = polar2xy(position_B)
    distance = np.sqrt(np.dot(positionA_xy-positionB_xy,positionA_xy-positionB_xy))
    return distance
def spiral(a,b,k):
    theta = (k*np.pi)/6
    position = np.array([theta,a*np.exp(b*theta)])
    return position
def angle_calculator(position_A, position_B,position_C):
    result = np.square(polar_distance(position_A, position_B))-np.square(polar_distance(position_C, position_A))-np.square(polar_distance(position_C, position_B))
    result = result / (-2*polar_distance(position_C, position_A)*polar_distance(position_B, position_C))
    result = np.arccos(result)
    return (result/np.pi)*180

# the log spiral is r = a*exp(b*theta)
a = 20
b = 0.1766
k= -22 # k is a integer, the theta is k*pi/6, where the circle is divided into 12 sectors.

print (block(a,b,k))
print (block(a,b,k-1))
print (block(a,b,k-2))
print (block(a,b,k-3))



print(angle_calculator(spiral(a,b,k+13),spiral(a,b,k+0),spiral(a,b,k+1)))