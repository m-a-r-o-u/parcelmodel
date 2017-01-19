import numpy as np
from state import State
from model import Model

def create_uniform_particle_distribution(total, groups):
    particle_count = (total / groups, ) * groups
    r_min = np.sort(np.random.uniform(1.e-6, 2.e-6, groups))
    return r_min, particle_count

def create_explicit_particle_distribution(r_min, particle_count):
    return r_min, particle_count

PARTICLE_DISTRIBUTIONS = {
    'uniform': create_uniform_particle_distribution,
    'explicit': create_explicit_particle_distribution,
    }

def create_particle_distribution(definitions):
    if definitions['type'] == 'uniform':
        total = float(definitions['total'])
        groups = int(definitions['groups'])
        return PARTICLE_DISTRIBUTIONS[definitions['type']](total, groups)
    elif definitions['type'] == 'explicit':
        r_min = map(float, definitions['r_min'])
        particle_count = map(float, definitions['particle_count'])
        return PARTICLE_DISTRIBUTIONS[definitions['type']](r_min, particle_count)

def model_init(initial_conditions):
    r_min, particle_count = create_particle_distribution(initial_conditions['particle_distribution'])
    model_parameters = {
        'r_min': r_min,
        'particle_count': particle_count,
        'T': initial_conditions['environment']['T'],
        'dt': initial_conditions['dt'],
        't_max': initial_conditions['t_max'],
        }
    intial_state = State(0, 280., 100000., 0.01, (0,)*len(r_min))
    return Model(model_parameters, intial_state)

