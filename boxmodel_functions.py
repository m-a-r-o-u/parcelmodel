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

def cloud_water(N, r):
  '''Return could water [kg kg-1] from number density N [m-3] and radius r [m]'''
  qc = 4. / 3. * np.pi * r ** 3 * c.RHO_H2O / c.RHO_AIR * N
  return qc

def radius(qc, N):
  '''Return mean radius in [m] from qc [kg kg-1] and N [m-3]'''
  r = (3. / 4. / np.pi * qc * c.RHO_AIR / c.RHO_H2O / N) ** (1./3.)
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

def conservative_gauss_perturbations(std, number):
    perturbations = np.random.normal(0., std, number)
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

def thermal_radiative_cooling_rate_using_libRadTran():
    '''Thermal Cloud Top heating rates in [K/s] for a cloud between 1-2 [km] of LWC=0.3 [g/m3], 
    reff = 10 [mum] average cooling -10 [K/d], resolved in 5 [m] steps'''
    hr = np.array([3.00082986111e-06, 2.61034375e-06, 2.27114699074e-06, 1.97577430556e-06, 1.71971296296e-06, 1.49758564815e-06, 1.30437962963e-06, 1.13560046296e-06, 9.89227314815e-07, 8.61581018519e-07, 7.50488194444e-07, 6.54801041667e-07, 5.7063912037e-07, 4.96938541667e-07, 4.33622569444e-07, 3.78971990741e-07, 3.30041782407e-07, 2.87311805556e-07, 2.5109212963e-07, 2.19143518519e-07, 1.90941666667e-07, 1.67326967593e-07, 1.46006712963e-07, 1.27233564815e-07, 1.11271782407e-07, 9.72445486111e-08, 8.45648263889e-08, 7.37137847222e-08, 6.56687268519e-08, 5.661125e-08, 4.88103703704e-08, 4.32703935185e-08, 3.76636458333e-08, 3.34929398148e-08, 2.92131481481e-08, 2.51227199074e-08, 2.18931597222e-08, 1.93660416667e-08, 1.67767013889e-08, 1.44392592593e-08, 1.25981828704e-08, 1.15278854167e-08, 1.00693530093e-08, 8.35482407407e-09, 7.46319675926e-09, 6.95775925926e-09, 5.63357986111e-09, 4.63761111111e-09, 4.62784722222e-09, 3.75544097222e-09, 3.13130671296e-09, 2.78786342593e-09, 2.34137037037e-09, 1.93591782407e-09, 2.01459606481e-09, 1.69188078704e-09, 9.5846087963e-10, 9.73941666667e-10, 6.86375694444e-10, 4.132875e-10, 6.71205324074e-10, -3.01757060185e-10, -1.04758101852e-10, 6.47805208333e-11, -5.58949884259e-10, -5.47386342593e-10, -1.18668171296e-09, -1.63712268519e-09, -1.11463784722e-09, -1.63551851852e-09, -2.21330671296e-09, -2.63371875e-09, -2.94757291667e-09, -3.50260300926e-09, -4.10700231481e-09, -4.42807523148e-09, -5.24622569444e-09, -5.93607407407e-09, -6.87014351852e-09, -7.75582638889e-09, -9.01591087963e-09, -1.06250543981e-08, -1.21111458333e-08, -1.38219907407e-08, -1.54176157407e-08, -1.72261574074e-08, -2.06154861111e-08, -2.41124652778e-08, -2.72358796296e-08, -3.09467476852e-08, -3.48639930556e-08, -4.02148726852e-08, -4.68529050926e-08, -5.30060648148e-08, -6.06023842593e-08, -6.99000694444e-08, -8.02931365741e-08, -9.18907060185e-08, -1.05061770833e-07, -1.20028472222e-07, -1.38431134259e-07, -1.59155324074e-07, -1.82180902778e-07, -2.08566550926e-07, -2.39625694444e-07, -2.7486400463e-07, -3.15388310185e-07, -3.62246990741e-07, -4.15166087963e-07, -4.76631481481e-07, -5.47413194444e-07, -6.2885e-07, -7.22248148148e-07, -8.29260300926e-07, -9.52400462963e-07, -1.09441909722e-06, -1.2570462963e-06, -1.44540625e-06, -1.6613275463e-06, -1.90874537037e-06, -2.19433333333e-06, -2.5233125e-06, -2.90188773148e-06, -3.33773842593e-06, -3.83966550926e-06, -4.41790625e-06, -5.0840625e-06, -5.85108449074e-06, -6.73585532407e-06, -7.75500115741e-06, -8.92991203704e-06, -1.028571875e-05, -1.18487962963e-05, -1.36530092593e-05, -1.57346759259e-05, -1.81370486111e-05, -2.09111574074e-05, -2.41145833333e-05, -2.78160069444e-05, -3.20923842593e-05, -3.7035625e-05, -4.27520601852e-05, -4.93627893519e-05, -5.70121527778e-05, -6.58649537037e-05, -7.61143518519e-05, -8.79864814815e-05, -0.000101742800926, -0.000117689699074, -0.00013618287037, -0.000157640509259, -0.000182548842593, -0.000211478356481, -0.000245101157407, -0.000284201967593, -0.000329708912037, -0.000382726041667, -0.000444563425926, -0.000516805324074, -0.000601388310185, -0.000700729861111, -0.000817944907407, -0.000957201736111, -0.00112435185185, -0.0013280775463, -0.00158199074074, -0.00190856134259, -0.00234653935185, -0.0029652025463, -0.00167974074074])
    return hr[-2] # [K/s]
#    kelvin_per_day = 250.
#    return - kelvin_per_day / (60.*60.*24.)

def thermal_radiation_using_libRadTran():
    hr = thermal_radiative_cooling_rate_using_libRadTran()
    return - hr * c.C_P / ( np.pi / 2) / 1.e8 / 1.e-10

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

def condensation(T, p, qv, qc_sum, qc, particle_count, r_min, dt, radiation, S_perturbation, math=np):
    r_old = math.maximum(r_min, radius(qc, particle_count))
    es = saturation_pressure(T)
    S = relative_humidity(T, p, qv) - 1 + S_perturbation
    if radiation:
       #E = thermal_radiation(T, qc_sum)
        E = thermal_radiation_using_libRadTran()
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
