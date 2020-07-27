import numpy as np
from .markov_schema import Markov_schema


def conservative_gauss_turbulence(Nsd, std):
    perturbations = np.random.normal(0., std, Nsd)
    return perturbations - perturbations.mean()


def conservative_gauss_turbulence_wrapper(t):
    def _f(Nsd, std, dt):
        def _g(r, N):
            return t(Nsd, std)
        return _g
    return _f


def no_turbulence(Nsd):
    return np.zeros(Nsd)


def no_turbulence_wrapper(t):
    def _f(Nsd, l, dt):
        def _g(r, N):
            return t(Nsd)
        return _g
    return _f


TURBULENCE_SCHEMES = {
    'no_turbulence': no_turbulence_wrapper(no_turbulence),
    'gauss_turbulence': conservative_gauss_turbulence_wrapper(conservative_gauss_turbulence),
    'markov_schema': Markov_schema,
}


def choose_turbulence_schema(definitions):
    definitions_t = dict(definitions['turbulence_schema'])
    definitions_t.update(
        {k: v for k, v in definitions.items() if k in ['Nsd', 'dt', 'l']})
    kwargs = {k: v for k, v in definitions_t.items() if k != 'type'}
    return TURBULENCE_SCHEMES[definitions_t['type']](**kwargs)
