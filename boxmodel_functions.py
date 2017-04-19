import boxmodel_constants as c
import numpy as np
import sh
from tempfile import NamedTemporaryFile
from scipy.integrate import odeint
from scipy.interpolate import interp1d

def saturation_pressure(T, math=np):
    '''Return saturation pressure [Pa] over flat water surface from temperature [K] valid only between 228.15 - 333.15'''
    return c.ES0 * math.exp(17.62 * (T - c.T0) / (243.12 + (T - c.T0)))

def saturation_vapor(T, p, math=np):
    '''Return saturation mixing ratio [kg kg-1] from T [K] and p [Pa]'''
    es = saturation_pressure(T, math=math)
    return c.R_G / c.R_V * es / (p - es)

def relative_humidity(T,p,qv):
    return qv / saturation_vapor(T,p)

def _cloud_water(N, r, rho=c.RHO_AIR):
    return  4. / 3. * np.pi * r ** 3 * c.RHO_H2O / rho * N

def cloud_water(N, r, r_min=0., rho=c.RHO_AIR):
    '''could water mixing ratio [kg kg-1] from number density N [# m-3] and radius r [m]'''
    return _cloud_water(N, r, rho=rho) - _cloud_water(N, r_min, rho=rho)

def _radius(qc, N, rho=c.RHO_AIR):
    return (3. / 4. / np.pi * qc * rho / c.RHO_H2O / N) ** (1./3.)

def radius(qc ,N ,r_min=0, rho=c.RHO_AIR):
    '''Return mean radius in [m] from qc [kg kg-1] and N [m-3]'''
    return _radius(qc + cloud_water(N, r_min, rho=rho) , N, rho=rho)

def stefan_boltzmann_law(T):
    '''Return power of a black body [W m-2] from T [K]'''
    return c.SIGMA_SB * T ** 4

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

def effective_radius(qc, N, r_min):
    '''Returns a effective radius because of r_min even if qc=0'''
    return np.sum(radius(qc, N, r_min) ** 3) / np.sum(radius(qc ,N ,r_min) ** 2)

def dynamic_cooling(w):
    lapse_rate = 10
    return -lapse_rate / 1000. * w

def fall_speed(r):
    '''approximation found in rogers: Short Course in Cloud Physics p.126'''
    assert np.all(r >= 0)
    assert np.all(r < 2.e-3)
    k1 = 1.19e8
    k2 = 8e3
    k3 = 2.01e2
    m1 = r < 40.e-6
    m3 = r > 0.6e-3
    m2 = ~(m1 | m3)
    r1 = k1 * r**2
    r2 = k2 * r
    r3 = k3 * r**0.5
    return -(m1*r1 + m2*r2 + m3*r3)

def differential_growth_by_condensation(r, t, E_net, es, T, S):
    '''differential diffusional growth equation returning dr/dt [m s-1] from ...units...'''
    c1 = c.H_LAT ** 2 / (c.R_V * c.K * T ** 2) + c.R_V * T / (c.D * es)
    c2 = c.H_LAT / (c.R_V * c.K * T ** 2)
    return (S / r + c2 * E_net) / (c1 * c.RHO_H2O)

def differential_growth_by_condensation_jacobian(r, t, E_net, es, T, S):
    '''jacobian of differential diffusional growth equation returning dr/dt [m s-1] from ...units...'''
    c1 = c.H_LAT ** 2 / (c.R_V * c.K * T ** 2) + c.R_V * T / (c.D * es)
    c2 = c.H_LAT / (c.R_V * c.K * T ** 2)
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
    return delta_T, delta_qv, delta_qc

def condensation_solver_euler(r_old, dt, E, es, T, S):
    return r_old + dt * differential_growth_by_condensation(r_old, dt, E, es, T, S)

def interp_afglus(file_name):
    z, p, T, air, o3, o2, h2o, co2, no2 = np.loadtxt(file_name, unpack=True)
    z *= 1000
    p *= 100
    def _interp_afglus(_z):
        pressure = interp1d(z, p)
        temperature = interp1d(z, T)
        return {'z':_z, 'p':pressure(_z), 'T':temperature(_z)}
    return _interp_afglus

def cloud_quantities(z, dz, qc, N, r_min):
    zmin = np.min(z) - np.min(z)%dz
    zmax = np.max(z) + (dz - np.max(z)%dz)
    zgrid = np.arange(zmin, zmax+dz, dz)
    weights = qc
    weights2 = radius(qc, N, r_min) ** 2
    weights3 = radius(qc, N, r_min) ** 3
    qc, bin_edge = np.histogram(z, zgrid, weights=weights)
    r2, bin_edge = np.histogram(z, zgrid, weights=weights2)
    r3, bin_edge = np.histogram(z, zgrid, weights=weights3)
    m = qc > 0
    r3[m] = r3[m] / r2[m]
    qc = np.append(qc, 0.0)
    r3 = np.append(r3, 0.0)
    return zgrid, qc, r3

def libRadTran_radiation_wrapper_thermal(height_in, lwc_in, reff_in):
    height = height_in[::-1] / float(1000)
    lwc = lwc_in[::-1] * float(1000)
    reff = reff_in[::-1] * 1.e6

    cloud_input = '\n'.join(['{}\t{}\t{}'.format(h, l, r) for h, l, r in zip(height, lwc, reff)])
    with NamedTemporaryFile(delete=False) as cloud_file:
        cloud_file.write(cloud_input)
        cloud_file.flush()
        cloud_file_name = cloud_file.name
        #what is layer_fd doing?
        uvspec_input_param = {
            'data_files_path': '/home/m/Mares.Barekzai/Software/libRadtran/data ',
            'atmosphere_file': '/home/m/Mares.Barekzai/Software/libRadtran/data/atmmod/afglus.dat',
            'source': 'thermal',
            'source_file': '/home/m/Mares.Barekzai/Software/libRadtran/data/solar_flux/fu',
            'albedo': '0',
            'rte_solver': 'disort',
            'number_of_streams': '4',
            'wavelength_index': '7 18',
            'mol_abs_param': 'FU',
            'wc_file_dimension': '1D',
            'wc_file': cloud_file_name,
            'output_process': 'sum',
            'zout': ' '.join((height_in/float(1000)).astype(str)),
            'layer___': 'layer_fd',
            }
        uvspec_input = """
        data_files_path {data_files_path}
        atmosphere_file {atmosphere_file}
        source {source} {source_file}
        albedo {albedo}
        rte_solver {rte_solver}
        number_of_streams {number_of_streams}
        wavelength_index {wavelength_index}
        mol_abs_param {mol_abs_param}
        wc_file {wc_file_dimension} {wc_file}
        output_process {output_process}
        zout {zout}
        heating_rate {layer___}
        quiet
        """.format(**uvspec_input_param)

        UVSPEC = sh.Command('uvspec')
        try:
            res = UVSPEC(_in=uvspec_input).stdout
        except sh.ErrorReturnCode as e:
            print "UVSPEC INPUT WAS:"
            print uvspec_input
            print "UVSPEC STDERR WAS:"
            print e.stderr
            print "UVSPEC STDOUT WAS:"
            print e.stdout
            raise
        return np.array(map(float, res.split())).reshape(len(height_in),2).T

def link_hight_to_radiation(z, zgrid, r, hr, rho_sp, rho=1):
    m = []
    for i in range(0, len(zgrid) - 1, 1):
        mlow = z < zgrid[i]
        mhig = z < zgrid[i+1]
        m.append(~mlow & mhig)
    res = np.zeros(len(z))
    for i in range(len(m)):
        res[m[i]] = heating_rate_to_Enet(hr[i], r[m[i]], rho_sp, rho=rho)
    return res

def heating_rate_to_Enet(hr, r, rho_sp, rho=1):
    r2_total = np.sum(r**2)
    return - hr * c.C_P * rho / 24. / 60. / 60. / r2_total / np.pi / rho_sp

def radiation_using_uvspec(z, dz, qc, particle_count, r_min):
    zgrid, lwc, reff = cloud_quantities(z, dz, qc, particle_count, r_min)
    res = libRadTran_radiation_wrapper_thermal(zgrid, lwc, reff)
    return link_hight_to_radiation(z, zgrid, radius(qc, particle_count, r_min), res[1], particle_count[0])
