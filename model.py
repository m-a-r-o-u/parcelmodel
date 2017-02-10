import numpy as np
import boxmodel_functions as bf
import boxmodel_constants as bc
from state import State

class Model(object):
    units = {
        't':'s',
        'T':'K',
        'p':'Pa',
        'qv':'kg kg-1',
        'qc':'kg kg-1',
    }

    def __init__(self, model_parameters, initial_state):
        self.r_min = np.array(model_parameters['r_min'])
        self.particle_count = np.array(model_parameters['particle_count'])
        self.radiation = model_parameters.get('radiation', False)
        self.T_env = model_parameters['T']
        self.dt = model_parameters['dt']
        self.t_max = model_parameters['t_max']
        initial_state = initial_state.copy()
        initial_state.qc = np.array(initial_state.qc, dtype='float')
        self._initial_state = initial_state
        self.distribution = model_parameters['distribution']
        self.distribution['radiation'] = self.radiation
        assert len(self.r_min) == len(self.particle_count)
        assert len(self.r_min) == len(self._initial_state.qc)

    def run(self, logger):
        state = self.initial_state()
        logger.set_units(self.units)
        logger.set_informations(self.distribution)
        logger.log_state(state)
        while not self.is_converged(state):
            state = self.step(state)
            logger.log_state(state)

    def initial_state(self):
        return self._initial_state.copy()

    def is_converged(self, state):
        return state.t >= self.t_max

    def step(self, old_state):
        new_state = self.prepare_new_state(old_state)
        delta_Ts, delta_qvs, delta_qc = self.calculate_tendencies(new_state)

        new_state.qc = new_state.qc + delta_qc
        new_state.T += np.sum(delta_Ts)
        new_state.qv += np.sum(delta_qvs)
        return new_state

    def calculate_tendencies(self, state):
        qc_sum = sum(state.qc)
        def condensation(qc, particle_count, r_min):
            return bf.condensation(state.T, state.p, state.qv, qc_sum, qc, particle_count, r_min, self.dt, self.radiation)
        delta_Ts, delta_qvs, delta_qc = condensation(state.qc, self.particle_count, self.r_min)
        return delta_Ts, delta_qvs, delta_qc

    def prepare_new_state(self, old_state):
        new_state = old_state.copy()
        new_state.t += self.dt
        qc_sum = sum(old_state.qc)
        cooling_rate = bf.thermal_radiative_cooling_rate(old_state.T, qc_sum, self.T_env)
        new_state.T += cooling_rate * self.dt
        return new_state
