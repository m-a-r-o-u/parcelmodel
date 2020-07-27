import numpy as np
from ccn import double_log_normal

ratio = 0.6
mu1 = 1.e-8
sigma1 = 1.4
mu2 = 7.5e-8
sigma2 = 1.6
Nsd = 1000

r, N  = double_log_normal(ratio, mu1, sigma1, mu2, sigma2, 1.e8, Nsd)
print np.max(r)
