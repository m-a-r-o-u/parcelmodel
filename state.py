class State(object):
    def __init__(self, t, T, p, qv, qc, z):
        self.t = t
        self.T = T
        self.p = p
        self.qv = qv
        self.qc = qc
        self.z = z
    
    def __repr__(self):
        return 'State({}, {}, {}, {}, {})'.format(self.t, self.T, self.p, self.qv, self.qc, self.z)

    def copy(self):
        return State(self.t, self.T, self.p, self.qv, self.qc, self.z)
