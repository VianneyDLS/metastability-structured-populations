import unittest
import numpy as np
from simulate import esd, energy, growth, competition, stochastic, deterministic

class ModelTests(unittest.TestCase):
    def test_esd(self):
        for A in (.5,1.5):
            x,m=esd(A)
            np.testing.assert_allclose(growth(x,A)-competition(x[:,None],x[None,:])@m,0,atol=1e-12)
            self.assertAlmostEqual(energy(A,x,m),0,places=12)
            z=np.linspace(-1,1,1001)
            self.assertLessEqual(np.max(growth(z,A)-competition(z[:,None],x[None,:])@m),1e-12)

    def test_stochastic(self):
        t=np.arange(11,dtype=float); K=80
        for A in (.5,1.5):
            x,n,E,mass,div,events=stochastic(A,K,t,t,123)
            self.assertEqual(n[0].sum(),K//2)
            self.assertTrue(np.all(n>=0))
            self.assertTrue(np.all(np.diff(div)<=0))
            self.assertFalse(np.any((n[:-1]==0)&(n[1:]>0)))
            np.testing.assert_allclose(n.sum(axis=1)/K,mass,atol=1e-12)
            np.testing.assert_allclose([energy(A,x,row/K) for row in n],E,atol=1e-11)
            again=stochastic(A,K,t,t,123)
            np.testing.assert_array_equal(n,again[1])

    def test_deterministic(self):
        for A in (.5,1.5):
            x,f,E=deterministic(A,np.arange(21,dtype=float),101)
            np.testing.assert_allclose(f[0],.25)
            self.assertTrue(np.all(np.diff(E)<=1e-9))
            self.assertGreaterEqual(f.min(),-1e-9)
            self.assertGreaterEqual(E.min(),-1e-9)

if __name__=='__main__':
    unittest.main()
