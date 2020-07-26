import numpy as np
import boxmodel_constants as c
import sh
from tempfile import NamedTemporaryFile
from boxmodel_functions import radius


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
            print("UVSPEC INPUT WAS:")
            print(uvspec_input)
            print("UVSPEC STDERR WAS:")
            print(e.stderr)
            print("UVSPEC STDOUT WAS:")
            print(e.stdout)
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

def thermal_radiation_using_uvspec(dz, state, microphysics):
    r_min = microphysics['r_min']
    particle_count = microphysics['particle_count']
    z = state.z
    qc = state.qc
    zgrid, lwc, reff = cloud_quantities(z, dz, qc, particle_count, r_min)
    res = libRadTran_radiation_wrapper_thermal(zgrid, lwc, reff)
    return link_hight_to_radiation(z, zgrid, radius(qc, particle_count, r_min), res[1], particle_count[0])
