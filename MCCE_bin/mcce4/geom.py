#!/usr/bin/env python
# geometry functions
import sys
import numpy as np
import math

# Geometry operations are through a 4x4 operation matrix.
# Operations can be stacked in the matrix.
# The matrix then is applied to a series of points.

class LINE:
    def __init__(self):
        self.p0 = np.array([0, 0, 0])   # line pass this  point
        self.t = np.array([1, 0, 0])    # direction cosines

    def from2p(self, p0, p1):
        self.p0 = np.array(p0)
        # self.p1 = np.array(p1)
        v = np.subtract(p1, p0)
        v_norm = np.linalg.norm(v)
        self.t = v/v_norm

    def print_me(self):
        print("Line passes (%.3f, %.3f, %.3f) in direction = (%.3f, %.3f, %.3f)" %
              (self.p0[0], self.p0[1], self.p0[2], self.t[0], self.t[1], self.t[2]))


class PLANE:
    def __init__(self) -> None:
        self.p0 = np.array([0, 0, 0])  # plane passes this point
        self.t = np.array([1, 0, 0])   # direction cosine of the plane (normal line)

    def from_3p(self, v1, v2, v3):
        """
        Return a plane object from 3 points.
        This plane passes 3 given points.
        The normal direction is defined by the cross product of direction v1->v2 and the direction 
        of v2 -> v3 
        """ 
        self.p0 = v1
        v12 = np.subtract(v2, v1)
        v23 = np.subtract(v3, v2)
        self.t = vector_normalize(np.cross(v12, v23)) 

class OPERATION:
    def __init__(self):
        self.operation = np.zeros((4, 4))
        self.operation[0][0] = 1.0
        self.operation[1][1] = 1.0
        self.operation[2][2] = 1.0
        self.operation[3][3] = 1.0
        return

    def reset(self):
        self.operation.fill(0)
        self.operation[0][0] = 1.0
        self.operation[1][1] = 1.0
        self.operation[2][2] = 1.0
        self.operation[3][3] = 1.0
        return

    def printme(self):
        print(self.operation)
        return

    def move(self, v):
        "Move by vector"
        op = np.zeros((4, 4))
        op[0][0] = 1.0
        op[1][1] = 1.0
        op[2][2] = 1.0
        op[3][3] = 1.0
        op[0][3] = v[0]
        op[1][3] = v[1]
        op[2][3] = v[2]
        self.operation = np.matmul(op, self.operation)

    def clone(self):
        clone_op = OPERATION()
        clone_op.operation = self.operation.copy()
        return clone_op

    def inverse(self):
        op = OPERATION()
        op.operation = np.linalg.inv(self.operation)
        return op

    def transpose(self):
        self.operation = self.operation.transpose()

    def roll(self, phi, axis):
        # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

        # validate the axis
        v = np.copy(axis.t)
        if np.linalg.norm(v) > 0.0000001:  # validate direction cosines
            v = v / np.linalg.norm(v)

            # translate to origin
            rotate_op = np.zeros((4, 4))
            rotate_op[0][0] = 1.0
            rotate_op[1][1] = 1.0
            rotate_op[2][2] = 1.0
            rotate_op[3][3] = 1.0
            rotate_op[0][3] = -axis.p0[0]
            rotate_op[1][3] = -axis.p0[1]
            rotate_op[2][3] = -axis.p0[2]
            self.operation = np.matmul(rotate_op, self.operation)

            # rotate
            SIN = np.sin(phi)
            COS = np.cos(phi)
            rotate_op = np.zeros((4, 4))
            rotate_op[3][3] = 1.0

            rotate_op[0][0] = v[0]*v[0]*(1-COS) + COS
            rotate_op[0][1] = v[0]*v[1]*(1-COS) - v[2]*SIN
            rotate_op[0][2] = v[0]*v[2]*(1-COS) + v[1]*SIN

            rotate_op[1][0] = v[0]*v[1]*(1-COS) + v[2]*SIN
            rotate_op[1][1] = v[1]*v[1]*(1-COS) + COS
            rotate_op[1][2] = v[1]*v[2]*(1-COS) - v[0]*SIN

            rotate_op[2][0] = v[0]*v[2]*(1-COS) - v[1]*SIN
            rotate_op[2][1] = v[1]*v[2]*(1-COS) + v[0]*SIN
            rotate_op[2][2] = v[2]*v[2]*(1-COS) + COS

            self.operation = np.matmul(rotate_op, self.operation)

            # translate back
            rotate_op = np.zeros((4, 4))
            rotate_op[0][0] = 1.0
            rotate_op[1][1] = 1.0
            rotate_op[2][2] = 1.0
            rotate_op[3][3] = 1.0
            rotate_op[0][3] = axis.p0[0]
            rotate_op[1][3] = axis.p0[1]
            rotate_op[2][3] = axis.p0[2]
            self.operation = np.matmul(rotate_op, self.operation)

        return

    def apply(self, v):  # replaces geom_apply()
        v1 = np.append(np.array(v), 1)[np.newaxis].transpose()
        v2 = np.matmul(self.operation, v1)
        v3 = v2[:-1]
        return v3.T[0]


# def geom_apply(operation, v):
#     v1 = np.append(np.array(v), 1)[np.newaxis].transpose()
#     v2 = np.matmul(operation.operation, v1)
#     v3 = v2[:-1]
#     return v3.T[0]

def vector_sum3v(v1, v2, v3):
    return((v1[0] + v2[0] + v3[0], v1[1] + v2[1] + v3[1], v1[2] + v2[2] + v3[2]))

def vector_scale(v, s):
    return((v[0]*s, v[1]*s, v[2]*s))

def vector_vminusv(v1, v2):
    return((v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]))

def vector_vplusv(v1, v2):
    return((v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]))

def vector_vxv(v1, v2):
    return((v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]))

def vector_normalize(v):
    #print(v)
    return v/np.linalg.norm(np.array(v))

def vector_orthogonal(v):
    if v[0] == 0 and v[1] == 0:
        return (0, 1, 0)
    else:
        return vector_normalize((-v[1], v[0], 0))

def ddvv(xyz1, xyz2):
    """Distance squared between two vectors."""
    dx=xyz1[0]-xyz2[0]
    dy=xyz1[1]-xyz2[1]
    dz=xyz1[2]-xyz2[2]
    return dx*dx+dy*dy+dz*dz

def dvv(xyz1, xyz2):
    """Distance  between two vectors."""
    return math.sqrt(ddvv(xyz1, xyz2))


def avv(v1, v2):
    "Angle between v1 to v2 in [0, pi], around axis with right-hand rule."
    v1_n = np.linalg.norm(v1)
    v2_n = np.linalg.norm(v2)
    COS = np.dot(v1, v2) / v1_n / v2_n

    return np.arccos(COS)

def geom_2v_onto_2v(v1, v2, t1, t2):
    "superimpose v1 to t1, then align v2 to t2"
    op = OPERATION()

    # superimpose v1 to t1
    op.move(np.array(t1) - np.array(v1))

    # align v1-v2 to t1-t2
    line1 = LINE()
    line1.from2p(v1, v2)
    line2 = LINE()
    line2.from2p(t1, t2)
    angle = avv(line1.t, line2.t)
    axis = LINE()
    axis.p0 = np.array(t1)
    axis.t = np.cross(line1.t, line2.t)
    op.roll(angle, axis)
    return op

def geom_3v_onto_3v(v1, v2, v3, t1, t2, t3):
    "superimpose v1 to t1, then align v2 to t2, lastly rotate to put v3 on the same plane as t3"
    # step 1 and 2, superimpose v1 to t1, align v2, to t2
    op = geom_2v_onto_2v(v1, v2, t1, t2)

    # step 3, align the normal direction of plane v1-v2-v3 to t1-t2-t3 
    # update
    new_v1 = op.apply(v1)
    new_v2 = op.apply(v2)
    new_v3 = op.apply(v3)
    plane_v123 = PLANE()
    plane_v123.from_3p(new_v1, new_v2, new_v3)
    plane_t123 = PLANE()
    plane_t123.from_3p(t1, t2, t3)
    # align
    angle = avv(plane_v123.t, plane_t123.t)
    axis = LINE()
    axis.p0 = t1
    
    # new_t = vector_vminusv(t2, t1)  # use t1->t2 as rotational axis direction
    # axis.t = vector_normalize(new_t)
    # Use t1->t2 has a problem to determine the direction of the angle.

    # Didn't handle 0 and 180 case yet
    axis.t = vector_normalize(np.cross(plane_v123.t, plane_t123.t))
    op.roll(angle, axis)

    return op



if __name__ == "__main__":
    p1 = (1, 1, 1)
    p2 = (9, 3, 3)

    t1 = (5, 1, 0)
    t2 = (0, 2, 0)

    op = geom_2v_onto_2v(p1, p2, t1, t2)
    op2 = op.inverse()

    new_p2 = op.apply(p2)
    recovered = op2.apply(new_p2)
    print("Test geom_2v_onto_2v:")
    print("Align %s, %s to %s, %s" % (str(p1), str(p2), str(t1), str(t2)))
    print("New p2:", new_p2)
    print("Recovered:", recovered)
    print()

    p1 = (2.364, -0.26, 0.041)
    p2 = (1.138, 0.515, 0.453)
    p3 = (3.01, 0.096, -0.916)

    t1 = (-10.322, 6.657, -13.101)
    t2 = (-9.007, 7.347, -13.224)
    t3 = (-10.316, 5.503, -12.577)
    op = geom_3v_onto_3v(p1, p2, p3, t1, t2, t3)
    new_p1 = op.apply(p1)
    new_p2 = op.apply(p2)
    new_p3 = op.apply(p3)

    print("Test geom_3v_onto_3v:")
    print("Align p %s, %s, %s to t %s, %s, %s" % (str(p1), str(p2), str(p3), str(t1), str(t2), str(t3)))
    print("New p1, p2, p3:", new_p1, new_p2, new_p3)
    op2 = op.inverse()
    recovered_p1 = op2.apply(new_p1)
    recovered_p2 = op2.apply(new_p2)
    recovered_p3 = op2.apply(new_p3)
    print("Recovered:", recovered_p1, recovered_p2, recovered_p3)
    print()
