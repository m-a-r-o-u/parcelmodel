import numpy as np
from state import State
from model import Model

def create_uniform_particle_distribution(total, groups):
    total = float(total)
    groups = int(groups)
    particle_count = (total / groups, ) * groups
    r_min = np.sort(np.random.uniform(1.e-6, 2.e-6, groups))
    return r_min, particle_count

def create_explicit_particle_distribution(r_min, particle_count):
    return map(float, r_min), map(float, particle_count)

def create_marine_ccn_particle_distribution(file_name, total, groups):
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

PARTICLE_DISTRIBUTIONS = {
    'uniform': create_uniform_particle_distribution,
    'explicit': create_explicit_particle_distribution,
    'marine_ccn': create_marine_ccn_particle_distribution,
    }

def create_particle_distribution(definitions):
    kwargs = {k:v for k,v in definitions.iteritems() if k != 'type'}
    return PARTICLE_DISTRIBUTIONS[definitions['type']](**kwargs)

def model_init(initial_conditions, executer):
    r_min, particle_count = create_particle_distribution(initial_conditions['particle_distribution'])
    model_parameters = {
        'r_min': r_min,
        'particle_count': particle_count,
        'T': initial_conditions['environment']['T'],
        'dt': initial_conditions['dt'],
        't_max': initial_conditions['t_max'],
        'radiation': initial_conditions['radiation'],
        'distribution': initial_conditions['particle_distribution'],
        'output_step': initial_conditions['output_step'],
        'perturbation': initial_conditions['perturbation'],
        'std': initial_conditions['std']
        }
    intial_state = State(0,
                         initial_conditions['T'],
                         initial_conditions['p'],
                         initial_conditions['qv'],
                         (0,)*len(r_min))
    return Model(model_parameters, intial_state, executer)
