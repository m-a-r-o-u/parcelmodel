'''hold frequently used constants of cloud microphysics'''


class Constant(float):
    '''Class for pysical constants'''

    def __new__(cls, value, units, doc):
        self = float.__new__(cls, value)
        self.units = units
        self.doc = doc
        return self


D            = Constant( 2.82e-5,     'm2 s-1',         'diffusion constant of water vapor in air')
K            = Constant( 2.43e-2,  'W m-1 K-1',                      'thermal conductivity of air')
R_V          = Constant( 461.401, 'J kg-1 K-1',             'specific gas constant of water vapor')
R_G          = Constant( 287.102, 'J kg-1 K-1',                               'ideal gas constant')
C_P          = Constant(  1003.5, 'J kg-1 K-1', 'specific heat of air at constant pressure')
C_P_H2O      = Constant(  4181.3, 'J kg-1 K-1', 'specific heat of water at constant pressure')
RHO_H2O      = Constant(   1000.,     'kg m-3',                                 'density of water')
RHO_AIR      = Constant(      1.,     'kg m-3',                                   'density of air')
H_LAT        = Constant( 2.257e6,     'J kg-1',                             'heat of vaporization')
SIGMA_SB     = Constant( 5.67e-8,  'W m-2 K-4',                        'stefan boltzmann constant')
ES0          = Constant(   611.2, 'kg m-1 s-2','saturation pressure over flat water at 273.15 [K]')
T0           = Constant(  273.15,          'K',                    'freezing temperature of water')
GAMMA        = Constant( 72.7e-3,      'N m-1',              'surface tension of water at 293 [K]')
N_A          = Constant(6.022e23,      'mol-1',           'particles per mol: avogradros constant')
M_MOL_H2O    = Constant(  18.e-3,   'kg mol-1',                          'molecular mass of water')
M_MOL_S      = Constant( 58.4e-3,   'kg mol-1',                           'molecular mass of NaCl')
M_MOL_AIR    = Constant(28.96e-3,   'kg mol-1',                    'average molecular mass of air')
RHO_S        = Constant(  2.16e3,     'kg m-3',                       'density of NaCL at 298 [K]')
ETA_AIR      = Constant( 17.1e-6,       'Pa s',              'dynamic viscosity of air at 273 [K]')
G            = Constant(    9.81,  'kg m2 s-2',                                 'gravity constant')
LAPSE_RATE_A = Constant( G / C_P,      'K m-1',                             'adiabatic lapse rate')
H_PLANCK     = Constant(6.62e-34,        'J s',                     'Plancksches Wirkungsquantuum')
K_BOLTZ      = Constant( 1.3e-23,      'J / K',                              'Boltzmann-Konstante')
C            = Constant(2.9979e8,      'm / s',                                   'Speed of Light')
S0           = Constant(    1367,      'W/ m2', 'solar constant')
