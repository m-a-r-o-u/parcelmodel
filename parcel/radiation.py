import numpy as np
from .boxmodel_functions import stefan_boltzmann_law
from .radiation_libradtran import thermal_radiation_using_uvspec


def stefan_boltzmann_schema(state, microphysics, factor, dz):
    E_net = stefan_boltzmann_law(state.T) * factor
    return E_net
    # return np.array([E_net] * len(state.qc))


def stefan_boltzmann_wrapper(f):
    def _f(factor, l):
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
        else:
            print('libradtran_wrapper mode Error')
    return _f


RADIATION_SCHEMES = {
    'libradtran': libradtran_wrapper(thermal_radiation_using_uvspec),
    'stefan_boltzmann': stefan_boltzmann_wrapper(stefan_boltzmann_schema),
    'no_radiation': no_radiation_wrapper(no_radiation),
}


def choose_radiation_schema(definitions):
    definitions_r = dict(definitions['radiation_schema'])
    definitions_r.update({k: v for k, v in definitions.items() if k in ['l']})
    kwargs = {k: v for k, v in definitions_r.items() if k != 'type'}
    return RADIATION_SCHEMES[definitions_r['type']](**kwargs)
