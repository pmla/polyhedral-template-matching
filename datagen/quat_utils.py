import numpy as np


def flip_up(q):
    if q[0] < 0:
        q = -q
    return q


def multiply(r, a):

    b0 = r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]
    b1 = r[0] * a[1] + r[1] * a[0] + r[2] * a[3] - r[3] * a[2]
    b2 = r[0] * a[2] - r[1] * a[3] + r[2] * a[0] + r[3] * a[1]
    b3 = r[0] * a[3] + r[1] * a[2] - r[2] * a[1] + r[3] * a[0]
    return flip_up(np.array([b0, b1, b2, b3]))


def random_quaternions(n):

    x0 = np.random.uniform(0, 1, n)
    x1 = np.random.uniform(0, 2 * np.pi, n)
    x2 = np.random.uniform(0, 2 * np.pi, n)

    r1, r2 = np.sqrt(1 - x0), np.sqrt(x0)
    s1, c1 = np.sin(x1), np.cos(x1)
    s2, c2 = np.sin(x2), np.cos(x2)

    return np.array([s1 * r1, c1 * r1, s2 * r2, c2 * r2]).T


def multiply_first_component(r, a):
    return r[0] * a[0] - r[1] * a[1] - r[2] * a[2] - r[3] * a[3]


def rotate_into_fundamental_zone(q, generators):

    index = np.argmax([abs(multiply_first_component(q, g)) for g in generators])
    return multiply(q, generators[index])


def map_points_out(basis_points, basis_weights, superset, subset, map_indices):

    assert( len(map_indices) * len(subset) == len(superset) )

    superset = np.array(superset)
    subset = np.array(subset)

    mapped_points = []
    mapped_weights = []
    for g in superset[map_indices]:
        for b, w in zip(basis_points, basis_weights):
            r = multiply(b, g)
            #r = rotate_into_fundamental_zone(r, subset)
            mapped_points += [r]
            mapped_weights += [w]

    return np.array(mapped_points), np.array(mapped_weights)


def quaternion_to_rotation_matrix(q):

    a, b, c, d = q

    u0 = a*a + b*b - c*c - d*d
    u1 = 2*b*c - 2*a*d
    u2 = 2*b*d + 2*a*c

    u3 = 2*b*c + 2*a*d
    u4 = a*a - b*b + c*c - d*d
    u5 = 2*c*d - 2*a*b

    u6 = 2*b*d - 2*a*c
    u7 = 2*c*d + 2*a*b
    u8 = a*a - b*b - c*c + d*d

    return np.array([[u0, u1, u2], [u3, u4, u5], [u6, u7, u8]])


def rodrigues_to_quaternion(r):

	s = 1 / np.sqrt(1 + np.linalg.norm(r)**2)
	return np.array([s, s * r[0], s * r[1], s * r[2]])


def rotation_matrix_to_quaternion(u):

    r11, r12, r13 = u[0]
    r21, r22, r23 = u[1]
    r31, r32, r33 = u[2]

    q0 = (1.0 + r11 + r22 + r33) / 4.0
    q1 = (1.0 + r11 - r22 - r33) / 4.0
    q2 = (1.0 - r11 + r22 - r33) / 4.0
    q3 = (1.0 - r11 - r22 + r33) / 4.0
    q = np.array([q0, q1, q2, q3])
    q = np.sqrt(np.maximum(0, q))

    i = np.argmax(q)
    if i == 0:
        q[1] *= np.sign(r32 - r23)
        q[2] *= np.sign(r13 - r31)
        q[3] *= np.sign(r21 - r12)

    elif i == 1:
        q[0] *= np.sign(r32 - r23)
        q[2] *= np.sign(r21 + r12)
        q[3] *= np.sign(r13 + r31)

    elif i == 2:
        q[0] *= np.sign(r13 - r31)
        q[1] *= np.sign(r21 + r12)
        q[3] *= np.sign(r32 + r23)

    elif i == 3:
        q[0] *= np.sign(r21 - r12)
        q[1] *= np.sign(r31 + r13)
        q[2] *= np.sign(r32 + r23)

    return q / np.linalg.norm(q)


#euler angle conventions are taken from EMSoft
def quaternion_to_euler(q):

    qq = q**2
    q03 = qq[0] + qq[3]
    q12 = qq[1] + qq[2]
    chi = np.sqrt(q03 * q12)
    if chi == 0:
        if q12 == 0:
            phi = 0
            phi2 = 0
            phi1 = np.arctan2(-2 * q[0] * q[3], qq[0] - qq[3])
        else:
            phi = np.pi
            phi2 = 0
            phi1 = np.arctan2(2 * q[1] * q[2], qq[1] - qq[2])
    else:
        phi = np.arctan2(2 * chi, q03 - q12)
        phi1 = np.arctan2((-q[0] * q[2] + q[1] * q[3]) / chi, (-q[0] * q[1] - q[2] * q[3]) / chi)
        phi2 = np.arctan2((+q[0] * q[2] + q[1] * q[3]) / chi, (-q[0] * q[1] + q[2] * q[3]) / chi)

    res = np.array([phi1, phi, phi2]) % (2 * np.pi)
    res[1] %= np.pi
    return res


#euler angle conventions are taken from EMSoft
def euler_to_quaternion(e):

    ee = np.array(e) / 2
    cphi = np.cos(ee[1])
    sphi = np.sin(ee[1])
    cm = np.cos(ee[0] - ee[2])
    sm = np.sin(ee[0] - ee[2])
    cp = np.cos(ee[0] + ee[2])
    sp = np.sin(ee[0] + ee[2])

    res = np.array([cphi * cp, -sphi * cm, -sphi * sm, -cphi * sp])
    return flip_up(res)

