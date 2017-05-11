import numpy as np
from boxmodel_functions import effective_radius
from boxmodel_functions import stefan_boltzmann_law
from boxmodel_functions import radius
from radiation_libradtran import thermal_radiation_using_uvspec

def optical_thickness(qc, microphysics):
    '''Stephens 1978b'''
    return 3. / 2. * np.sum(qc) / effective_radius(qc, microphysics['particle_count'], microphysics['r_min'])

def flux_at_drop(flux, n, r, dz):
    return flux / np.pi / dz / np.sum(n * r**2)

def thermal_radiation(state, microphysics, math=np):
    tau = optical_thickness(state.qc, microphysics)
    return stefan_boltzmann_law(state.T) * (1 - math.exp(-tau))

def stefan_boltzmann_schema(state, microphysics, factor, dz):
    E_net = np.array([thermal_radiation(state, microphysics)] * len(state.qc))
    r = radius(state.qc, microphysics['particle_count'], microphysics['r_min'])
    n = microphysics['particle_count']
    return np.minimum(flux_at_drop(E_net, n, r, dz), E_net) * factor

def stefan_boltzmann_wrapper(f):
    def _f(factor=1, l=100):
        def stefan_boltzmann_schema(state, microphysics):
            return f(state, microphysics, factor, l)
        return stefan_boltzmann_schema
    return _f

def no_radiation(state, microphysics):
    return np.zeros(len(state.qc))

def no_radiation_wrapper(t):
    def _f():
        def no_radiation_schema(state, microphysics):
            return t(state, microphysics)
        return no_radiation_schema
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
    'stefan_boltzmann': stefan_boltzmann_wrapper(stefan_boltzmann_schema),
    'no_radiation': no_radiation_wrapper(no_radiation),
    }

def choose_radiation_schema(definitions):
    definitions_r = dict(definitions['radiation_schema'])
    definitions_r.update({k:v for k, v in definitions.iteritems() if k in ['l'] })
    kwargs = {k:v for k,v in definitions_r.iteritems() if k != 'type'}
    return RADIATION_SCHEMES[definitions_r['type']](**kwargs)
