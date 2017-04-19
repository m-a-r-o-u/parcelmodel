import numpy as np

class Markov_schema(object):
    def __init__(self, groups, l, epsilon, del_t):
        self.groups = groups
        self.l = l
        self.epsilon = epsilon
        self.w = np.zeros(groups)
        self.S = np.zeros(groups)
        self.del_t = del_t

    def __call__(self, r=np.array([10.e-6]), N=np.array([100.e6])):
        S = self.S
        self.S = self.saturation_fluctuations(r, N)
        tke = self.sgs_tke()
        self.w = self.ornstein_uhlenbeck_process(self.turbulent_timescale(tke), self.std_deviation_w(tke))
        return S - S.mean()
    
    def sgs_tke(self, Ce=0.845):
        return (self.l * self.epsilon/ Ce ) ** (2. / 3.)
    
    def turbulent_timescale(self, tke, Ct=1.5):
        return self.l / (2. * np.pi) ** (1. / 3.) * (Ct / tke) ** (1. / 2.)
    
    def std_deviation_w(self, tke):
        return (2. / 3. * tke) ** (1. / 2.)
    
    def ornstein_uhlenbeck_process(self, tau, sigma):
        return self.w * np.exp(-self.del_t / tau) + (1 - np.exp(-2. * self.del_t / tau)) ** (1. / 2.) * sigma * np.random.normal(0., 1., self.groups)
    
    def saturation_fluctuations(self, r, N):
        A1 = 3.e-4
        A2 = 2.8e-4
        tau_relax = (A2 * np.sum(r * N)) ** -1
        if tau_relax < self.del_t: raise Error
        return self.del_t * (A1 * self.w - self.S / tau_relax)
