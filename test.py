import unittest
from gaussian import GaussianInteger
from dgt import dgt_cooley_tukey, idgt_gentleman_sande
import random
from params import N, p

NRUNS = 10**2

def gen_polynomial_modp(length):
    x = []

    for _ in range(length):
        x.append(random.randrange(0,p))
    
    return x

def schoolbook_mul(a, b):
    assert len(a) == len(b)
    
    N = len(a)
    c = [0]*N

    for i in range(N):
        for j in range(N):
            v = a[i] * b[j] * (-1)**(int((i+j)//float(N)))
            c[(i + j) % N] = (c[(i + j) % N] + v) % p
    return c 

class TestNegacyclicDGT(unittest.TestCase):

    def test_negacyclic(self):
        print("\nDGT transformation merged with the right-angle convolution")

        for _ in range(NRUNS):            
            a = gen_polynomial_modp(N)            
            A = dgt_cooley_tukey(a)
            A = idgt_gentleman_sande(A)

            self.assertEqual(a, A)

    def test_mul_negacyclic(self):
        print("\nPolynomial multiplication using DGT merged with the right-angle convolution")

        for _ in range(NRUNS):
            a = gen_polynomial_modp(N)
            b = gen_polynomial_modp(N)

            ab = schoolbook_mul(a, b)

            a_dgt = dgt_cooley_tukey(a)
            b_dgt = dgt_cooley_tukey(b)

            c_dgt = [x * y for x, y in zip(a_dgt, b_dgt)]

            c = idgt_gentleman_sande(c_dgt)

            self.assertEqual(c, ab)    

if __name__ == '__main__':
    unittest.main()