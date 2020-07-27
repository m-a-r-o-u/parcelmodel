import numpy as np
from state import State
from model import Model
from ccn import choose_particle_distribution
from radiation import choose_radiation_schema
from turbulence import choose_turbulence_schema
from atmosphere import choose_atmosphere_schema


def model_init(initial_conditions, executer):
    r_min, particle_count = choose_particle_distribution(initial_conditions)
    radiation_function = choose_radiation_schema(initial_conditions)
    turbulence_schema = choose_turbulence_schema(initial_conditions)
    atmosphere_schema = choose_atmosphere_schema(initial_conditions)

    model_parameters = {
        'r_min': r_min,
        'particle_count': particle_count,
        'w': initial_conditions['w'],
        'l': float(initial_conditions['l']),
        'feedback': initial_conditions['feedback'],
        'dt': initial_conditions['dt'],
        't_max': initial_conditions['t_max'],
        'dt_output': max(initial_conditions['dt_output'], initial_conditions['dt']),
        'radiation_function': radiation_function,
        'turbulence_schema': turbulence_schema,
        'atmosphere_schema': atmosphere_schema,
        'distribution': initial_conditions['particle_distribution'],
    }
    initial_state = State(0,
                          initial_conditions['T'],
                          initial_conditions['p'],
                          initial_conditions['qv'],
                          np.array([0.] * len(r_min)),
                          np.array([initial_conditions['z0']] * len(r_min)),
                          0,
                          np.array([0.] * len(r_min)),
                          np.array([0.] * len(r_min))
                          )
    return Model(model_parameters, initial_state, executer)
