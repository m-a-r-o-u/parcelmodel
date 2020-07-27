import numpy as np

class Markov_schema(object):
    '''
    The Markov_schema object parametrizes saturation fluctuations

    :param Nsd: Initialize internal arrays, that hold S, w for each particle
    :type Nsd: Integer?
    :ivar Nsd: This is where we store Nsd
    '''
    def __init__(self, Nsd, l, epsilon, dt, gen = np.random.RandomState()):
        self.__name__ = self.__class__.__name__
        self.Nsd = Nsd
        self.l = l
        self.epsilon = epsilon
        self.w = np.zeros(Nsd)
        self.S = np.zeros(Nsd)
        self.dt = min(dt, 0.02)
        self.n_average = max(int(dt / 0.02), 1)
        self.gen = gen

    def __call__(self, r=np.array([10.e-6]), N=np.array([100.e6])):
        S_average = np.zeros(self.Nsd)
        for t in range(self.n_average):
             S_average += self.S_step(r=r, N=N) / self.n_average
        return S_average - S_average.mean()

    def S_step(self, r=np.array([10.e-6]), N=np.array([100.e6])):
        S = self.S
        self.S = S + self.saturation_fluctuations(r, N)
        tke = self.sgs_tke()
        self.w = self.ornstein_uhlenbeck_process(self.turbulent_timescale(tke), 
                                                 self.std_deviation_w(tke))
        return S
    
    def sgs_tke(self, Ce=0.845):
        return (self.l * self.epsilon/ Ce ) ** (2. / 3.)
    
    def turbulent_timescale(self, tke, Ct=1.5):
        return self.l / (2. * np.pi) ** (1. / 3.) * (Ct / tke) ** (1. / 2.)
    
    def std_deviation_w(self, tke):
        return (2. / 3. * tke) ** (1. / 2.)
    
    def ornstein_uhlenbeck_process(self, tau, sigma):
        return self.w * np.exp(-self.dt / tau) + (1 - np.exp(-2. * self.dt / tau)) ** (1. / 2.) * sigma * self.gen.normal(0., 1., self.Nsd)
    
    def saturation_fluctuations(self, r, N):
        A1 = 3.e-4
        A2 = 2.8e-4
        tau_relax = (A2 * np.sum(r * N)) ** -1
        if tau_relax < self.dt: raise Error
        return self.dt * (A1 * self.w - self.S / tau_relax)
