import numpy as np
from state import State
from model import Model
from ccn import choose_particle_distribution
from radiation import choose_radiation_schema
from turbulence import choose_turbulence_schema

def model_init(initial_conditions, executer):
    r_min, particle_count = choose_particle_distribution(initial_conditions)
    radiation_function = choose_radiation_schema(initial_conditions['radiation_schema'])
    turbulence_schema = choose_turbulence_schema(initial_conditions)
    model_parameters = {
        'r_min': r_min,
        'particle_count': particle_count,
        'T': initial_conditions['environment']['T'],
        'w': initial_conditions['environment']['w'],
        'dt': initial_conditions['dt'],
        't_max': initial_conditions['t_max'],
        'radiation_function': radiation_function,
        'turbulence_schema': turbulence_schema,
        'distribution': initial_conditions['particle_distribution'],
        'output_step': initial_conditions['output_step'],
        }
    initial_state = State(0,
                         initial_conditions['T'],
                         initial_conditions['p'],
                         initial_conditions['qv'],
                         np.array([0.] * len(r_min)),
                         np.array([initial_conditions['environment']['z0']] * len(r_min)),
                         np.array([0.] * len(r_min)))
    return Model(model_parameters, initial_state, executer)
