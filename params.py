from gaussian import GaussianInteger

p = 12289
N = 16
r = 11
k = N//2

g = pow(r, (p-1)//k, p)
assert(pow(g, k, p) == 1)
g_inv = pow(g, p-2, p)
assert((g * g_inv) % p == 1)

invkmodp = pow(k, p-2, p)

kth_root_of_i = GaussianInteger(3408, 7802)
kth_root_of_i_inv = GaussianInteger(2575, 2377)
assert(pow(kth_root_of_i, k) % p == GaussianInteger(0, 1))
assert(pow(kth_root_of_i_inv, k) % p == GaussianInteger(0, p-1))