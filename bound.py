import numpy as np
from scipy.special import binom as binom
from os import sys

m = int(sys.argv[1])
t = int(sys.argv[2])

cur_k = 2**(m-t)
cur_k_proj = 0
l = 1
k_lim = cur_k + binom(t, l) * 2**(m-t)
while l <= m:
    k_step = np.lcm(l, t) // l
    k_proj_step = np.lcm(l, t) // t
    while cur_k < k_lim:
        print(cur_k, end = " ")
        print(cur_k_proj)
        cur_k += k_step
        cur_k_proj += k_proj_step
    l += 1
    k_lim += binom(t, l) * 2**(m-t)