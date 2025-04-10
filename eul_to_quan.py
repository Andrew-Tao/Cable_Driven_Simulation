import numpy as np


def euler_to_quaternion(roll, pitch, yaw):
    """
    Convert Euler angles to a quaternion.
    :param roll: Rotation around x-axis in radians
    :param pitch: Rotation around y-axis in radians
    :param yaw: Rotation around z-axis in radians
    :return: Quaternion (x, y, z, w)
    """
    cy = np.cos(yaw * 0.5)
    sy = np.sin(yaw * 0.5)
    cp = np.cos(pitch * 0.5)
    sp = np.sin(pitch * 0.5)
    cr = np.cos(roll * 0.5)
    sr = np.sin(roll * 0.5)

    w = cr * cp * cy + sr * sp * sy
    x = sr * cp * cy - cr * sp * sy
    y = cr * sp * cy + sr * cp * sy
    z = cr * cp * sy - sr * sp * cy

    return (x, y, z, w)





def quaternion_to_euler(x, y, z, w):
    """
    Convert a quaternion into Euler angles (roll, pitch, yaw)
    Roll is rotation around the x-axis, pitch is rotation around the y-axis,
    and yaw is rotation around the z-axis.

    :param x: x component of the quaternion
    :param y: y component of the quaternion
    :param z: z component of the quaternion
    :param w: w component (scalar part) of the quaternion
    :return: tuple of Euler angles (roll, pitch, yaw) in radians
    """
    # Roll (x-axis rotation)
    sinr_cosp = 2 * (w * x + y * z)
    cosr_cosp = 1 - 2 * (x * x + y * y)
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    # Pitch (y-axis rotation)
    sinp = 2 * (w * y - z * x)
    if np.abs(sinp) >= 1:
        # Use 90 degrees if out of range (result in a cliff effect)
        pitch = np.pi / 2 * np.sign(sinp)
    else:
        pitch = np.arcsin(sinp)

    # Yaw (z-axis rotation)
    siny_cosp = 2 * (w * z + x * y)
    cosy_cosp = 1 - 2 * (y * y + z * z)
    yaw = np.arctan2(siny_cosp, cosy_cosp)

    return roll, pitch, yaw


# Example usage:
# Define an example quaternion (x, y, z, w)
q =   (0.0, -0.7071067811865475, 0.0, 0.7071067811865476)# Example values; ideally, the quaternion should be normalized

# Normalize the quaternion
norm = np.sqrt(sum(component ** 2 for component in q))
q_normalized = tuple(component / norm for component in q)

roll, pitch, yaw = quaternion_to_euler(*q_normalized)
print("Euler angles (radians):")
print("Roll:  ", roll)
print("Pitch: ", pitch)
print("Yaw:   ", yaw)
