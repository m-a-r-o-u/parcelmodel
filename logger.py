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

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.finalize()

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

class NetCDFLogger(BaseLogger):
    def __init__(self):
        self.states = []

    def log_state(self, state):
        self.states.append(state)

    def finalize(self):
        from netCDF4 import Dataset
        with Dataset('./states.nc', mode='w', format='NETCDF3_64BIT') as file_handle:
            t_dim_nc = 'time'
            particle_dim_nc = 'super_particles'
            file_handle.createDimension(t_dim_nc, len(self.states))
            file_handle.createDimension(particle_dim_nc, len(self.states[0].qc))
            time_nc = file_handle.createVariable('time', 'i4', (t_dim_nc))
            T_nc = file_handle.createVariable('T', 'f4', (t_dim_nc))
            p_nc = file_handle.createVariable('p', 'f4', (t_dim_nc))
            qv_nc = file_handle.createVariable('qv', 'f4', (t_dim_nc))
            qc_nc = file_handle.createVariable('qc', 'f4', (t_dim_nc, particle_dim_nc))

            time_nc[:] = [state.t for state in self.states]
            T_nc[:] = [state.T for state in self.states]
            p_nc[:] = [state.p for state in self.states]
            qv_nc[:] = [state.qv for state in self.states]
            qc_nc[:] = [state.qc for state in self.states]

LOGGERS = {
    'MultiPlotLogger': PlotTimeSeriesLogger,
    'Logger': Logger,
    'NetCDFLogger': NetCDFLogger,
     }
