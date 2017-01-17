class State(object):
    def __init__(self, T, p, qv, qc, iteration_count=0):
        self.T = T
        self.p = p
        self.qv = qv
        self.qc = qc
        self.iteration_count = iteration_count
    
    def __repr__(self):
        return 'State({}, {}, {}, {})'.format(self.T, self.p, self.qv, self.qc)

    def copy(self):
        return State(self.T, self.p, self.qv, self.qc, self.iteration_count)
