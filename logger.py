class Logger(object):
    def finalize(self):
        pass

    def log_state(self, state):
        print state

class FinalStateLogger(object):
    def __init__(self):
        self.last_state = None

    def log_state(self, state):
        self.last_state = state

    def finalize(self):
        print self.last_state

class PlotTLogger(object):
    def finalize(self):
        pass
    
    def log_state(self, state):
        print state.T
 
class MultiLogger(object):
    def __init__(self, loggers):
        self.loggers = loggers

    def log_state(self, state):
        for logger in self.loggers:
            logger.log_state(state)

    def finalize(self):
        for logger in self.loggers:
            logger.finalize()
