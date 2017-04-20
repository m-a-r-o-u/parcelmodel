import numpy as np
from boxmodel_functions import effective_radius
from boxmodel_functions import stefan_boltzmann_law
from radiation_libradtran import thermal_radiation_using_uvspec

def optical_thickness(qc, microphysics):
    '''Stephens 1978b'''
    return 3. / 2. * np.sum(qc) / effective_radius(qc, microphysics['particle_count'], microphysics['r_min'])

def thermal_radiation(state, microphysics, T_env=250, math=np):
    tau = optical_thickness(state.qc, microphysics)
    return (stefan_boltzmann_law(state.T) - stefan_boltzmann_law(T_env)) * (1 - math.exp(-tau))

#???CHECK???
def thermal_radiative_cooling_rate(T, qc, N, r_min, T_env=250.):
    '''cooling rate due to radiation [units???] from cloud water mixing ratio [kg kg-1]'''
    E_net = thermal_radiation(T, qc, N, r_min, T_env=T_env)
    cooling_rate = - E_net * 80 / 20 / (3600. * 24.)
    return cooling_rate

def stefan_boltzmann_schema(state, microphysics):
    #muss noch auf die drops verteilt werden
    E_net = np.array([thermal_radiation(state, microphysics)] * len(state.qc))
    return E_net

def no_radiation(state, microphysics):
    return np.zeros(len(state.qc))

def function_wrapper(t):
    def _f():
        def _g(state, microphysics):
            return t(state, microphysics)
        return _g
    return _f

def libradtran_wrapper(t):
    def _f(dz, mode):
        if mode == 'thermal':
            def _g(state, microphysics):
                return t(dz, state, microphysics)
            return _g
        else: print 'libradtran_wrapper mode Error'
    return _f

RADIATION_SCHEMES = {
    'libradtran': libradtran_wrapper(thermal_radiation_using_uvspec),
    'stefan_boltzmann': function_wrapper(stefan_boltzmann_schema),
    'no_radiation': function_wrapper(no_radiation),
    }

def choose_radiation_schema(definitions):
    kwargs = {k:v for k,v in definitions.iteritems() if k != 'type'}
    return RADIATION_SCHEMES[definitions['type']](**kwargs)
