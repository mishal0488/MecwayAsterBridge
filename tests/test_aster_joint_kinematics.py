import unittest
import numpy as np
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "Scripts"))
from aster_joint_kinematics import from_vector,multiply,matrix,constraint

class JointTests(unittest.TestCase):
    def test_allowed_swing_blocked_twist(self):
        a=from_vector([.7,-.4,.9]);u=np.array([1.,0.,0.])
        for v in ([0,.8,0],[0,0,.8],[0,.6,.8]):
            b=multiply(a,from_vector(v))
            self.assertAlmostEqual(constraint(a,b,u)[0],0,places=13)
        self.assertGreater(abs(constraint(a,multiply(a,from_vector([.3,0,0])),u)[0]),.2)

    def test_objectivity_and_equal_opposite_couples(self):
        a=from_vector([.4,.7,-.2]);b=from_vector([-.3,.9,.5]);u=np.array([1.,0.,0.])
        q=from_vector([.8,-1.1,.6]);c,g=constraint(a,b,u)
        c2,g2=constraint(multiply(q,a),multiply(q,b),u)
        self.assertAlmostEqual(c,c2,places=13)
        np.testing.assert_allclose(g[:3]+g[3:],0,atol=1e-14)
        np.testing.assert_allclose(g2[:3],matrix(q)@g[:3],atol=1e-13)
        np.testing.assert_allclose(g2[3:],matrix(q)@g[3:],atol=1e-13)

    def test_spatial_jacobian_by_finite_differences(self):
        rng=np.random.default_rng(418)
        for _ in range(25):
            a=from_vector(rng.normal(size=3)*.6);b=from_vector(rng.normal(size=3)*.6)
            u=rng.normal(size=3);u/=np.linalg.norm(u)
            c,g=constraint(a,b,u);d=rng.normal(size=6);h=1e-6
            cp=constraint(multiply(from_vector(h*d[:3]),a),multiply(from_vector(h*d[3:]),b),u)[0]
            cm=constraint(multiply(from_vector(-h*d[:3]),a),multiply(from_vector(-h*d[3:]),b),u)[0]
            self.assertAlmostEqual((cp-cm)/(2*h),np.dot(g,d),places=8)

    def test_no_twist_after_compound_common_rotations(self):
        a=from_vector([0,0,0]);b=from_vector([0,.6,.8]);u=np.array([1.,0.,0.])
        for turn in ([.8,0,0],[0,-.7,.3],[.6,.4,-.9]):
            q=from_vector(turn);a=multiply(q,a);b=multiply(q,b)
            self.assertAlmostEqual(constraint(a,b,u)[0],0,places=13)

    def test_quaternion_sign_and_chart_limit(self):
        a=from_vector([.2,.1,0]);b=from_vector([.4,.6,.1]);u=np.array([1.,0.,0.])
        c,g=constraint(a,b,u);c2,g2=constraint(-a,b,u)
        self.assertAlmostEqual(c,c2,places=13);np.testing.assert_allclose(g,g2)
        with self.assertRaises(ValueError):constraint(from_vector([0,0,0]),from_vector([0,np.pi,0]),u)

if __name__=='__main__':unittest.main()
