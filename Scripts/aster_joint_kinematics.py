# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Objective, zero-relative-twist V/W joint kinematics (Code_Aster runtime).

Unit quaternions are scalar-first and map initial global frames to current
frames. u is the beam's initial global U axis. The holonomic constraint is
2*u.dot(vector(conj(qa)*qb)) == 0, on the branch short of a 180 degree swing.
It retains zero swing-twist torsion and two bending freedoms. No penalty
stiffness is used. aster_vw_joint supplies the multiplier constraints through
an outer nonlinear solve around native STAT_NON_LINE load increments.
"""
import math
import numpy as np


def multiply(a,b):
    return np.r_[a[0]*b[0]-np.dot(a[1:],b[1:]), a[0]*b[1:]+b[0]*a[1:]+np.cross(a[1:],b[1:])]


def conjugate(q):
    return q*np.array([1.,-1.,-1.,-1.])


def from_vector(v):
    v=np.asarray(v,dtype=float);t=np.linalg.norm(v)
    scale=.5-t*t/48 if t<1e-7 else math.sin(t/2)/t
    return np.r_[math.cos(t/2),scale*v]


def matrix(q):
    w=q[0];v=q[1:]
    K=np.array([[0.,-v[2],v[1]],[v[2],0.,-v[0]],[-v[1],v[0],0.]])
    return np.eye(3)+2*w*K+2*K@K


def constraint(qa,qb,u):
    qa=np.asarray(qa,dtype=float);qb=np.asarray(qb,dtype=float)
    u=np.asarray(u,dtype=float)
    if abs(np.linalg.norm(qa)-1)>1e-10 or abs(np.linalg.norm(qb)-1)>1e-10:
        raise ValueError('Unit quaternions required')
    if abs(np.linalg.norm(u)-1)>1e-10:raise ValueError('Unit U axis required')
    r=multiply(conjugate(qa),qb)
    if abs(r[0])<1e-8:raise ValueError('180-degree relative rotation is outside this joint chart')
    sign=1. if r[0]>0 else -1.
    r*=sign
    # Spatial virtual rotations: qb <- exp(db)*qb, qa <- exp(da)*qa.
    # The equal/opposite moments ensure no net joint couple under a common
    # virtual rigid rotation. Differentiate C=2*u.vector(r).
    g=matrix(qa) @ (r[0]*u+np.cross(r[1:],u))
    return 2*np.dot(u,r[1:]),np.r_[-g,g]
