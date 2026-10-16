"""Counterexamples to extending fixed linear rotation-vector release equations.
This is a kinematic check, not a solver implementation of selective hinges.
"""
import itertools,json,math
import numpy as np

def exp(v):
    t=np.linalg.norm(v)
    if t<1e-14:return np.eye(3)
    a=v/t;K=np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
    return np.eye(3)+math.sin(t)*K+(1-math.cos(t))*(K@K)
def log(R):
    t=math.acos(np.clip((np.trace(R)-1)/2,-1,1))
    if t<1e-14:return np.zeros(3)
    return t/(2*math.sin(t))*np.array([R[2,1]-R[1,2],R[0,2]-R[2,0],R[1,0]-R[0,1]])
report={}
for size in (1,2):
    for free in itertools.combinations(range(3),size):
        retained=next(i for i in range(3) if i not in free)
        f=free[0];other=next(i for i in range(3) if i not in (retained,f))
        Q=exp(np.eye(3)[other]*1.1)
        relative=np.eye(3)[f]*.6
        B=Q@exp(relative)
        np.testing.assert_allclose(log(Q.T@B),relative,atol=1e-12)
        residual=log(B)-log(Q)
        assert abs(residual[retained])>.1
        report[''.join('uvw'[i] for i in free)]={'allowed_relative_rotation':relative.tolist(),'spurious_retained_component_radians':float(residual[retained])}
print(json.dumps(report,indent=2))
