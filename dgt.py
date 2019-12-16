from gaussian import GaussianInteger
from params import N, p, g, g_inv, invkmodp, kth_root_of_i, kth_root_of_i_inv
from math import log

fold   = lambda a: [GaussianInteger(x, y) for x, y in zip(a[:N//2], a[N//2:])]
unfold = lambda a: [a[j].re for j in range(N//2)] + [a[j].imag for j in range(N//2)]
bit_reverse = lambda n, width: int('{:0{width}b}'.format(n, width=width)[::-1], 2)
log_2 = lambda n: n.bit_length()

def dgt_cooley_tukey(X):

    n = len(X)
    X = fold(X)

    d = n//4
    m = 1
    while(m < n//2):
        for k in range(m):
            j1 = 2*k*d
            j2 = j1 + d - 1
            a = (pow(g, bit_reverse(k, log_2(n//4)-1), p) * pow(kth_root_of_i, d)) % p
            for j in range(j1, j2+1, 1):
                u = a * X[j+d]
                X[j+d] = (X[j] - u) % p
                X[j] = (X[j] + u) % p
        d = d//2
        m = 2*m

    return X

def idgt_gentleman_sande(X):

    n = len(X)*2

    d = 1
    m = n//2
    while(m > 1):
        for j1 in range(d):
            r = 0
            i = j1
            while(i < n//2):
                u = X[i]
                v = X[i+d]
                X[i] = (u + v) % p
                X[i+d] = (pow(g_inv, bit_reverse(r, log_2(n//4)-1), p) * pow(kth_root_of_i_inv, d) * (u - v)) % p
                r += 1
                i = i + 2*d
        d = d*2
        m = m//2

    return unfold([(xi * invkmodp) % p for i, xi in enumerate(X)])