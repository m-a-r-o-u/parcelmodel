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
        'z':'m',
        'E':'W m-2',
    }

    def __init__(self, model_parameters, initial_state, executer):
        self.r_min = np.array(model_parameters['r_min'])
        self.particle_count = np.array(model_parameters['particle_count'])
        self.radiation = model_parameters.get('radiation', False)
        self.w = model_parameters['w']
        self.T_env = model_parameters['T']
        self.dt = model_parameters['dt']
        self.dz = model_parameters['dz']
        self.t_max = model_parameters['t_max']
        initial_state = initial_state.copy()
        initial_state.qc = np.array(initial_state.qc, dtype='float')
        self._initial_state = initial_state
        self.output_step = model_parameters['output_step']
        self.perturbation = model_parameters.get('perturbation', False)
        self.std = model_parameters.get('std')
        self.information = model_parameters['distribution']
        self.information['radiation'] = self.radiation
        self.information['perturbation'] = self.perturbation
        self.information['std'] = self.std
        self.atmosphere = bf.interp_afglus('./input/afglus.dat')
        assert len(self.r_min) == len(self.particle_count)
        assert len(self.r_min) == len(self._initial_state.qc)
        self.executer = executer(step, self, initial_state, self.output_step)
        self.step = self.executer.step

    def run(self, logger):
        state = self.initial_state()
        logger.set_units(self.units)
        logger.set_informations(self.information)
        logger.log_state(state)
        while not self.is_converged(state):
            state = self.step(state)
            logger.log_state(state)

    def initial_state(self):
        return self._initial_state.copy()

    def is_converged(self, state):
        return state.t >= self.t_max

    def calculate_tendencies(self, state, math=np):
        qc_sum = math.sum(state.qc)
        S_perturbations = bf.conservative_gauss_perturbations(self.std, len(self.r_min), self.perturbation)
        m = nucleation_slice(state, S_perturbations, self.r_min, self.particle_count)
        def condensation(qc, particle_count, r_min, S_perturbation, E):
            return bf.condensation(state.T, state.p, state.qv, qc_sum, qc, particle_count, r_min, self.dt, E, S_perturbation, math=math)
        delta_Ts, delta_qvs, delta_qc = np.zeros((3, len(state.qc)))
        delta_Ts[m], delta_qvs[m], delta_qc[m] = condensation(state.qc[m],
                                                              self.particle_count[m],
                                                              self.r_min[m],
                                                              S_perturbations[m],
                                                              state.E[m])
        return delta_Ts, delta_qvs, delta_qc

    def prepare_new_state(self, old_state, math=np):
        new_state = old_state.copy()
        new_state.t += self.dt
        print new_state.t
        qc_sum = math.sum(old_state.qc)
        cooling_rate = bf.dynamic_cooling(self.w)
        new_state.z = new_state.z + self.w * self.dt
        new_state.p = self.atmosphere(new_state.z.mean())['p']
        new_state.T = self.atmosphere(new_state.z.mean())['T']
        new_state.E = bf.radiation_using_uvspec(new_state.z, self.dz, new_state.qc, self.particle_count, self.r_min)
        return new_state

def nucleation_slice(state, S_perturbations, r_min, particle_count):
    S = bf.relative_humidity(state.T, state.p, state.qv) - 1 + S_perturbations
    m = (S < bf.critical_super_saturation(r_min, state.T)) & (np.isclose(bf.radius(state.qc, particle_count, r_min), r_min))
    return [not i for i in m]

def step(model, old_state, math=np):
    new_state = model.prepare_new_state(old_state, math=math)
    delta_Ts, delta_qvs, delta_qc = model.calculate_tendencies(new_state, math=math)
    new_state.qc = new_state.qc + delta_qc
    new_state.T += math.sum(delta_Ts)
    new_state.qv += math.sum(delta_qvs)
    new_state.z = new_state.z + bf.fall_speed(bf.radius(new_state.qc, model.particle_count, model.r_min)) * model.dt
    return new_state
