import numpy as np
from boxmodel_functions import effective_radius
from boxmodel_functions import stefan_boltzmann_law

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
    return 0, E_net

def no_radiation(state, microphysics):
    return 0, np.zeros(len(state.qc))

RADIATION_SCHEMES = {
    'stefan_boltzmann': stefan_boltzmann_schema,
    'no_radiation': no_radiation,
    }

def choose_radiation_schema(definitions):
    return RADIATION_SCHEMES[definitions['type']]
