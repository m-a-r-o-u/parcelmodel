def logger_factory(config):
    if isinstance(config, list):
        return MultiLogger([logger_factory(c) for c in config])
    else:
        constructor = LOGGERS[config['type']]
        kwargs = {k:v for k,v in config.iteritems() if k != 'type'}
        return constructor(**kwargs)

class BaseLogger(object):
    def set_units(self, units):
        pass

    def finalize(self):
        pass

class Logger(BaseLogger):
    def log_state(self, state):
        print state

class FinalStateLogger(BaseLogger):
    def __init__(self):
        self.last_state = None

    def log_state(self, state):
        self.last_state = state

    def finalize(self):
        print self.last_state

class PlotTLogger(BaseLogger):
    def __init__(self):
        self.t = []
        self.T = []
    
    def log_state(self, state):
        self.t.append(state.t)
        self.T.append(state.T)

    def finalize(self):
        import matplotlib.pyplot as plt
        plt.plot(self.t, self.T)
        plt.savefig('T.png')

class PlotQVLogger(BaseLogger):
    def __init__(self):
        self.t = []
        self.qv = []
    
    def log_state(self, state):
        print state.T
        self.t.append(state.t)
        self.qv.append(state.qv)

    def finalize(self):
        import matplotlib.pyplot as plt
        plt.plot(self.t, self.qv)
        plt.savefig('qv.png')

class MultiLogger(BaseLogger):
    def __init__(self, loggers):
        self.loggers = loggers

    def set_units(self, units):
        for logger in self.loggers:
            logger.set_units(units)

    def log_state(self, state):
        for logger in self.loggers:
            logger.log_state(state)

    def finalize(self):
        for logger in self.loggers:
            logger.finalize()

class PlotTimeSeriesLogger(BaseLogger):
    def __init__(self, quantities):
        self.units = []
        self.t = []
        self.quantities = quantities
        self.data = [[] for _ in quantities]

    def set_units(self, units):
        self.units = [units[quantity] for quantity in self.quantities]

    def log_state(self, state):
        self.t.append(state.t)
        for name, storage in zip(self.quantities, self.data):
            storage.append(getattr(state, name))

    def finalize(self):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(self.quantities))
        for name, storage, ax, unit in zip(self.quantities, self.data, axes, self.units):
            ax.plot(self.t, storage)
            ax.set_ylabel('{} [{}]'.format(name, unit))
        axes[-1].set_xlabel('time')
        fig.savefig('daten.png')

LOGGERS = {
    'MultiPlot': PlotTimeSeriesLogger,
    'Logger': Logger,
     }
