import numpy as np

def uniform_particle_distribution(total, groups):
    total = float(total)
    groups = int(groups)
    particle_count = (total / groups, ) * groups
    r_min = np.sort(np.random.uniform(1.e-6, 2.e-6, groups))
    return r_min, particle_count

def explicit_particle_distribution(r_min, particle_count):
    return map(float, r_min), map(float, particle_count)

def marine_ccn_particle_distribution(file_name, total, groups):
    from scipy import interpolate
    xdata, ydata = np.loadtxt(file_name, unpack=True)
    ysum = np.sum(ydata)
    ycumsum = np.cumsum(ydata) / ysum
    yinterp = interpolate. interp1d(ycumsum, xdata)
    r_min = np.random.uniform(np.min(ycumsum), np.max(ycumsum), groups)
    total = float(total)
    groups = int(groups)
    particle_count = (total / groups, ) * groups
    return r_min*1.e-6, particle_count

def double_log_normal(ratio, mu1, sigma1, mu2, sigma2, total, groups):
    n1 = int(ratio * groups)
    n2 = groups - n1
    r_ccn1 = np.random.lognormal(np.log(mu1), sigma1, n1)
    r_ccn2 = np.random.lognormal(np.log(mu2), sigma2, n2)
    r_min = np.concatenate((r_ccn1, r_ccn2))
    particle_count = (total / groups, ) * groups
    rmx=1.e-6
    print("try out maximal r_ccn of: " + str(rmx))
    r_min[r_min>rmx] = rmx
    return r_min, particle_count

PARTICLE_DISTRIBUTIONS = {
    'uniform': uniform_particle_distribution,
    'explicit': explicit_particle_distribution,
    'marine_ccn': marine_ccn_particle_distribution,
    'double_log_normal': double_log_normal,
    }

def choose_particle_distribution(definitions):
    definitions_p = dict(definitions['particle_distribution'])
    definitions_p.update({k:v for k, v in definitions.iteritems() if k == 'groups' })
    kwargs = {k:v for k,v in definitions_p.iteritems() if k != 'type'}
    return PARTICLE_DISTRIBUTIONS[definitions_p['type']](**kwargs)
