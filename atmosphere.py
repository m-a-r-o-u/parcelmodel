import numpy as np
from scipy.interpolate import interp1d
import boxmodel_constants as c

def interp_afglus_wrapper(atmosphere_file):
    z, p, T, air, o3, o2, h2o, co2, no2 = np.loadtxt(atmosphere_file, unpack=True)
    z *= 1000
    p *= 100
    def interp_afglus(_p, _T, z_in, dz):
        pressure = interp1d(z, p)
        temperature = interp1d(z, T)
        z_out = z_in + dz
        return [pressure(z_in), temperature(z_in), z_out]
    return interp_afglus

def linear_and_hydrostatic_wrapper(lapse_rate=c.LAPSE_RATE_A, rho_air=c.RHO_AIR):
    def linear_and_hydrostatic(p, T, z, dz, lapse_rate=lapse_rate, rho_air=rho_air):
        p_out = p - rho_air * c.G * dz
        T_out = T - lapse_rate * dz
        z_out = z + dz
        return [p_out, T_out, z_out]
    return linear_and_hydrostatic

ATMOSPHERE_SCHEMES = {
    'afglus': interp_afglus_wrapper,
    'linear_and_hydrostatic': linear_and_hydrostatic_wrapper,
    }

def choose_atmosphere_schema(definitions):
    definitions_t = definitions['atmosphere_schema']
    kwargs = {k:v for k,v in definitions_t.iteritems() if k != 'type'}
    return ATMOSPHERE_SCHEMES[definitions_t['type']](**kwargs)
