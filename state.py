class State(object):
    def __init__(self, t, T, p, qv, qc, z, E, age):
        self.t = t
        self.T = T
        self.p = p
        self.qv = qv
        self.qc = qc
        self.z = z
        self.E = E
        self.age = age
    
    def __repr__(self):
        return 'State({}, {}, {}, {}, {}, {}, {}, {})'.format(self.t, self.T, self.p, self.qv, self.qc, self.z, self.E, self.age)

    def copy(self):
        return State(self.t, self.T, self.p, self.qv, self.qc, self.z, self.E, self.age)
