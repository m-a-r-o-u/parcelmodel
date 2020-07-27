from . import boxmodel_constants as c
import numpy as np


def potential_temperature(Tin, p, inverse=False):
    e = c.R_G / c.C_P
    if inverse:
        e = -e
    return Tin * (p[0] / p) ** e


def stokes_falling_speed(r):
    '''return falling speed for water sphere in air'''
    return 2. / 9. * r * r * c.G * (c.RHO_H2O - c.RHO_AIR) / c.ETA_AIR


def saturation_pressure(T, math=np):
    '''Return saturation pressure [Pa] over flat water surface from temperature [K] valid only between 228.15 - 333.15'''
    return c.ES0 * math.exp(17.62 * (T - c.T0) / (243.12 + (T - c.T0)))


def saturation_vapor(T, p, math=np):
    '''Return saturation mixing ratio [kg kg-1] from T [K] and p [Pa]'''
    es = saturation_pressure(T, math=math)
    return c.R_G / c.R_V * es / (p - es)


def relative_humidity(T, p, qv):
    return qv / saturation_vapor(T, p)


def _cloud_water(N, r, rho=c.RHO_AIR):
    return 4. / 3. * np.pi * r ** 3 * c.RHO_H2O / rho * N


def cloud_water(N, r, r_min=0., rho=c.RHO_AIR):
    '''could water mixing ratio [kg kg-1] from number density N [# m-3] and radius r [m]'''
    return _cloud_water(N, r, rho=rho) - _cloud_water(N, r_min, rho=rho)


def _radius(qc, N, rho=c.RHO_AIR):
    return (3. / 4. / np.pi * qc * rho / c.RHO_H2O / N) ** (1. / 3.)


def radius(qc, N, r_min=0, rho=c.RHO_AIR):
    '''Return mean radius in [m] from qc [kg kg-1] and N [m-3]'''
    return _radius(qc + cloud_water(N, r_min, rho=rho), N, rho=rho)


def stefan_boltzmann_law(T):
    '''Return power of a black body [W m-2] from T [K]'''
    return c.SIGMA_SB * T ** 4


def kelvins_parameter(T=273.15):
    '''
    Informations are taken from the book:
    A short Course in Cloud Physics
    by
    Rogers and Yau
    page: 89
    '''
    return 2 * c.GAMMA / c.R_V / c.RHO_H2O / T


def kelvin_curvature_effect(r, T=273.15, math=np):
    A = kelvins_parameter(T=T)
    return math.exp(A / r)


def raoults_parameter(r_ccn):
    return 2 * r_ccn ** 3 * c.RHO_S * c.M_MOL_H2O / c.RHO_H2O / c.M_MOL_S


def raoult_mixture_effect(r, r_ccn):
    B = raoults_parameter(r_ccn)
    return (1 - B / r ** 3)


def critical_radius(r_ccn, T=273.15, math=np):
    return math.sqrt(3. * raoults_parameter(r_ccn) / kelvins_parameter(T=T))


def critical_super_saturation(r_ccn, T=273.15, math=np):
    return math.sqrt(4. * kelvins_parameter(T=T) ** 3 / 27. / raoults_parameter(r_ccn))


def koehler(r, r_ccn, T=273.15, math=np):
    return kelvin_curvature_effect(r, T=T, math=math) * raoult_mixture_effect(r, r_ccn)


def effective_radius(qc, N, r_min):
    '''Returns a effective radius because of r_min even if qc=0'''
    return np.sum(radius(qc, N, r_min) ** 3) / np.sum(radius(qc, N, r_min) ** 2)


def dynamic_cooling(w):
    lapse_rate = 10
    return -lapse_rate / 1000. * w


def fall_speed(r):
    '''approximation found in rogers: Short Course in Cloud Physics p.126'''
    assert np.all(r >= 0), 'r is {}'.format(r)
    assert np.all(r < 2.e-3)
    k1 = 1.19e8
    k2 = 8e3
    k3 = 2.01e2

    m1 = r < 40.e-6
    m2 = r >= 40.e-6
    m2 &= r < 0.6e-3
    m3 = r > 0.6e-3

    r1 = k1 * r**2
    r2 = k2 * r
    r3 = k3 * r**0.5

    return m1 * r1 + m2 * r2 + m3 * r3


def c2_differential_growth_by_condensation(T, es):
    return c.H_LAT / (c.R_V * c.K * T ** 2)


def c1_differential_growth_by_condensation(T, es):
    return c.H_LAT ** 2 / (c.R_V * c.K * T ** 2) + c.R_V * T / (c.D * es)


def differential_growth_by_condensation(r, t, E_net, es, T, S):
    '''differential diffusional growth equation returning dr/dt [m s-1] from ...units...'''
    c1 = c1_differential_growth_by_condensation(T, es)
    c2 = c2_differential_growth_by_condensation(T, es)
    return (S / r + c2 * E_net) / (c1 * c.RHO_H2O)


def differential_growth_by_condensation_jacobian(r, t, E_net, es, T, S):
    '''jacobian of differential diffusional growth equation returning dr/dt [m s-1] from ...units...'''
    c1 = c1_differential_growth_by_condensation(T, es)
    return - S / r ** 2 / (c1 * c.RHO_H2O)


def condensation(T, p, qv, qc_sum, qc, particle_count, r_min, dt, E, S_perturbation, math=np):
    r_old = math.maximum(r_min, radius(qc, particle_count, r_min))
    es = saturation_pressure(T)
    S = relative_humidity(T, p, qv) - 1 + S_perturbation
    r_new = condensation_solver_euler(r_old, dt, E, es, T, S)
    r_new = math.maximum(r_new, r_min)
    delta_qc = cloud_water(particle_count, r_new, r_min) - qc
    delta_T = delta_qc * c.H_LAT / c.C_P
    delta_qv = -delta_qc
    return np.array(delta_T), np.array(delta_qv), np.array(delta_qc)


def condensation_solver_euler(r_old, dt, E, es, T, S):
    return r_old + dt * differential_growth_by_condensation(r_old, None, E, es, T, S)


def planck(l, T):
    '''power of planck radiation per wavelength in W/(m2 m)'''
    c1 = 2. * np.pi * c.H_PLANCK * c.C * c.C
    c2 = c.H_PLANCK * c.C / c.K_BOLTZ
    return c1 / (l**5) / (np.exp(c2 / (l * T)) - 1)


def diffusion_solution_classic(t, S, T, es):
    A = c1_differential_growth_by_condensation(T, es)
    return np.sqrt(2. * S * t / A / c.RHO_H2O)
