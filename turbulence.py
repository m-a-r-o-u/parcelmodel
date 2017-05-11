import numpy as np
from markov_schema import Markov_schema

def conservative_gauss_turbulence(groups, std):
    perturbations = np.random.normal(0., std, groups)
    return perturbations - perturbations.mean()

def conservative_gauss_turbulence_wrapper(t):
    def _f(groups, std, dt):
        def _g(r, N):
            return t(groups, std)
        return _g
    return _f

def no_turbulence(groups):
    return np.zeros(groups)

def no_turbulence_wrapper(t):
    def _f(groups, l, dt):
        def _g(r, N):
            return t(groups)
        return _g
    return _f

TURBULENCE_SCHEMES = {
    'no_turbulence': no_turbulence_wrapper(no_turbulence),
    'gauss_turbulence': conservative_gauss_turbulence_wrapper(conservative_gauss_turbulence),
    'markov_schema': Markov_schema,
    }

def choose_turbulence_schema(definitions):
    definitions_t = dict(definitions['turbulence_schema'])
    definitions_t.update({k:v for k, v in definitions.iteritems() if k in ['groups', 'dt', 'l'] })
    kwargs = {k:v for k,v in definitions_t.iteritems() if k != 'type'}
    return TURBULENCE_SCHEMES[definitions_t['type']](**kwargs)
