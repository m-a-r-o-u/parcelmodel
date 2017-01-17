'''hold frequently used constants of cloud microphysics'''

class Constant(float):
  '''Class for pysical constants'''
  def __new__(cls, value, units, doc):
    self = float.__new__(cls, value)
    self.units = units
    self.doc = doc
    return self

D        = Constant(2.82e-5,     'm2 s-1',         'diffusion constant of water vapor in air')
K        = Constant(2.43e-2,  'W m-1 K-1',                      'thermal conductivity of air')
R_V      = Constant(461.401, 'J kg-1 K-1',             'specific gas constant of water vapor')
R_G      = Constant(287.102, 'J kg-1 K-1',                               'ideal gas constant')
C_P      = Constant( 1003.5, 'J kg-1 K-1',               'specific heat at constant pressure')
RHO_H2O  = Constant(  1000.,     'kg m-3',                                 'density of water')
RHO_AIR  = Constant(     1.,     'kg m-3',                                   'density of air')
H_LAT    = Constant(2.257e6,     'J kg-1',                             'heat of vaporization')
SIGMA_SB = Constant(5.67e-8,  'W m-2 K-4',                        'stefan boltzmann constant')
ES0      = Constant(  611.2, 'kg m-1 s-2','saturation pressure over flat water at 273.15 [K]')
T0       = Constant( 273.15,          'K',             'freezing temperature of water in [K]')


#R       = 287.102           #[J/(kg K)]
#R_V     = 461.401           #[J/(kg K)]
#CP      = 1003.5            #[J/(kg K)]
#D       = 0.282e-4          #[m2 s-1]    diff. constante (water air)
#RHO_AIR = 1.                #[kg/m3]
#RHO_H2O = 1000.             #[kg/m3]
#H_LAT   = 2257.e3           #[J kg-1]    heat of vaporization
#K       = 0.0243            #[W m-1 K-1] thermal conductivity
