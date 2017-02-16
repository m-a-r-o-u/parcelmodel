import boxmodel_constants as c
import numpy as np
from scipy.integrate import odeint


# add second doc string, to the right location
def saturation_pressure(T, math=np):
  '''Return saturation pressure [Pa] over flat water surface from temperature [K]'''
  #T_min = 228.15
  #T_max = 333.15
  #assert (T_min < T < T_max), "{0} is out range ({1}, {2}) [K] for magnus approximation".format(T, T_min, T_max)
  es = c.ES0 * math.exp(17.62 * (T - c.T0) / (243.12 + (T - c.T0)))
  return es

def saturation_vapor(T, p, math=np):
  '''Return saturation mixing ratio [kg kg-1] from T [K] and p [Pa]'''
  es = saturation_pressure(T, math=math)
  qvs = c.R_G / c.R_V * es / (p - es)
  return qvs

def saturation_temperature(p, qv, math=np):
  '''Return saturation temperature [K] from p [Pa] and qv[kg kg-1]'''
  x  = qv * c.R_V * p / c.R_G / c.ES0
  c1 = 1. / c.T0  
  c2 = c.R_V / c.H_LAT 
  T  = (c1 - math.log(x) * c2) ** (-1)
  return T

def saturation_adjustment(T, qv, qc, p, math=np):
  '''Return the state T [K], qv [kg kg-1], qv [kg kg-1] after saturation adjustment'''
  delqc = math.maximum(0, qv - bf.saturation_vapor(T,p))
  qc += delqc
  qv -= delqc
  T  += delqc * bc.H_LAT / bc.C_P
  return T,qv,qc

def relative_humidity(T,p,qv):
  '''return the relative humidity'''
  r = qv / saturation_vapor(T,p)
  return r

def cloud_water(N, r, math=np):
  '''Return could water [kg kg-1] from number density N [m-3] and radius r [m]'''
  qc = 4. / 3. * math.pi * r ** 3 * c.RHO_H2O / c.RHO_AIR * N
  return qc

def radius(qc, N, math=np):
  '''Return mean radius in [m] from qc [kg kg-1] and N [m-3]'''
  r = (3. / 4. / math.pi * qc * c.RHO_AIR / c.RHO_H2O / N) ** (1./3.)
  return r

def stefan_boltzmann_law(T):
  '''Return power of a black body [W m-2] from T [K]'''
  P = c.SIGMA_SB * T ** 4
  return P

def kelvins_parameter(T=273.15):
    return 2 * c.GAMMA / c.R_V / c.RHO_H2O / T

def kelvin_curvature_effect(r, T=273.15, math=np):
    A = kelvins_parameter(T=T)
    return math.exp(A / r )

def raoults_parameter(r_ccn):
    return 2 * r_ccn ** 3  * c.RHO_S * c.M_MOL_H2O / c.RHO_H2O / c.M_MOL_S

def raoult_mixture_effect(r, r_ccn):
    B = raoults_parameter(r_ccn)
    return (1 - B / r ** 3)

def critical_radius(r_ccn, T=273.15, math=np):
    return math.sqrt(3. * raoults_parameter(r_ccn) / kelvins_parameter(T=T))

def critical_super_saturation(r_ccn, T=273.15, math=np):
    return math.sqrt(4. * kelvins_parameter(T=T) ** 3 / 27. / raoults_parameter(r_ccn))

def koehler(r, r_ccn, T=273.15, math=np):
    return kelvin_curvature_effect(r, T=T, math=math) * raoult_mixture_effect(r, r_ccn)

def conservative_gauss_perturbations(std, number, math=np):
    perturbations = math.random.normal(0., std, number)
    return perturbations - perturbations.mean()

#???CHECK???
def optical_thickness(qc):
  '''optical thickness [units???] from cloud water mixing ratio [kg kg-1]'''
  tau = qc * 10. * 1000. + 0.1
  return tau

#???CHECK???
def thermal_radiation(T, qc, T_env=250, math=np):
  '''thermal radiation [units???] from cloud water mixing ratio [kg kg-1]'''
  tau   = optical_thickness(qc)
  E_net = (stefan_boltzmann_law(T) - stefan_boltzmann_law(T_env))* (1 - math.exp(-tau))
  return E_net

#???CHECK???
def thermal_radiative_cooling_rate(T, qc, T_env=250.):
  '''cooling rate due to radiation [units???] from cloud water mixing ratio [kg kg-1]'''
  E_net = thermal_radiation(T, qc, T_env=T_env)
  cooling_rate = - E_net * 80 / 20 / (3600. * 24.)
  return cooling_rate

def differential_growth_by_condensation(r, t, E_net, es, T, S):
  '''differential diffusional growth equation returning dr/dt [m s-1] from ...units...'''
  c1 = c.H_LAT ** 2 / (c.R_V * c.K * T ** 2) + c.R_V * T / (c.D * es)
  c2 = c.H_LAT / (c.R_V * c.K * T ** 2)
  r_new = (S / r + c2 * E_net) / (c1 * c.RHO_H2O)
  return r_new

def differential_growth_by_condensation_jacobian(r, t, E_net, es, T, S):
  '''jacobian of differential diffusional growth equation returning dr/dt [m s-1] from ...units...'''
  c1 = c.H_LAT ** 2 / (c.R_V * c.K * T ** 2) + c.R_V * T / (c.D * es)
  c2 = c.H_LAT / (c.R_V * c.K * T ** 2)
  r_new = - S / r ** 2 / (c1 * c.RHO_H2O)
  return r_new

def condensation(T, p, qv, qc_sum, qc, particle_count, r_min, dt, radiation, math=np):
    r_old = math.maximum(r_min, radius(qc, particle_count))
    es = saturation_pressure(T)
    S = relative_humidity(T, p, qv) - 1
    if radiation:
        E = thermal_radiation(T, qc_sum)
    else:
        E = 0
    r_new = condensation_solver_linear(r_old, dt, E, es, T, S)
    r_new = math.maximum(r_new, r_min)
    delta_qc = cloud_water(particle_count, r_new) - qc
    delta_T = delta_qc * c.H_LAT / c.C_P
    delta_qv = -delta_qc
    return delta_T, delta_qv, delta_qc

def condensation_solver_linear(r_old, dt, E, es, T, S):
    return r_old + dt * differential_growth_by_condensation(r_old, dt, E, es, T, S)
