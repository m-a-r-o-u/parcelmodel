class Logger(object):
    def finalize(self):
        pass

    def log_state(self, state):
        print state

class FinalStateLogger(object):
    def __init__(self):
        self.last_state = None

    def finalize(self):
        print self.last_state

    def log_state(self, state):
        self.last_state = state

class PlotTLogger(object):
    def finalize(self):
        pass
    
    def log_state(self, state):
        print state.T
