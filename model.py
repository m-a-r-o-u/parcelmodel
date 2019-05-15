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
        'age':'s',
    }

    def __init__(self, model_parameters, initial_state, executer):
        self.microphysics = {'r_min' : np.array(model_parameters['r_min']),
                             'particle_count' : np.array(model_parameters['particle_count'])}
        self.schemata = {'radiation': model_parameters['radiation_function'],
                         'turbulence': model_parameters['turbulence_schema'],
                         'atmosphere': model_parameters['atmosphere_schema']}
        self.w = model_parameters['w']
        self.l = model_parameters['l']
        self.feedback = model_parameters['feedback']
        self.dt = model_parameters['dt']
        self.t_max = int(model_parameters['t_max'])
        initial_state = initial_state.copy()
        initial_state.qc = np.array(initial_state.qc, dtype='float')
        self._initial_state = initial_state
        self.dt_output = model_parameters['dt_output']
        assert len(self.microphysics['r_min']) == len(self.microphysics['particle_count'])
        assert len(self.microphysics['r_min']) == len(self._initial_state.qc)
        self.executer = executer(step, self, initial_state, int(self.dt_output / float(self.dt)))
        self.step = self.executer.step

    def run(self, logger):
        state = self.initial_state()
        logger.set_units(self.units)
        logger.log_state(state)
        while not self.is_converged(state):
            print state.t
            state = self.step(state)
            logger.log_state(state)

    def initial_state(self):
        return self._initial_state.copy()

    def is_converged(self, state):
        return state.t >= self.t_max

    def calculate_tendencies(self, state, math=np):
        qc_sum = math.sum(state.qc)
        if qc_sum > 0.0:
            S_perturbations = self.schemata['turbulence'](bf.radius(state.qc,
                                                                self.microphysics['particle_count'],
                                                                self.microphysics['r_min']),
                                                      self.microphysics['particle_count'])
        else:
            S_perturbations = np.zeros(len(state.qc))


        def condensation(qc, particle_count, r_min, S_perturbation, E):
            return bf.condensation(state.T, state.p, state.qv, qc_sum, qc, particle_count, r_min, self.dt, E, S_perturbation, math=math)
        m = nucleation_slice(state, S_perturbations, self.microphysics)
        delta_Ts, delta_qvs, delta_qc = np.zeros((3, len(state.qc)))
        delta_Ts[m], delta_qvs[m], delta_qc[m] = condensation(state.qc[m],
                                                              self.microphysics['particle_count'][m],
                                                              self.microphysics['r_min'][m],
                                                              S_perturbations[m],
                                                              state.E)

        delta_Ts_E = -(state.E * self.dt / bc.C_P / bc.RHO_AIR / self.l)
        delta_age = np.array([0.]*len(self.microphysics['particle_count']))
        delta_age[m] += self.dt
        return delta_Ts * self.feedback['latent_heat'], delta_qvs, delta_qc, delta_Ts_E * self.feedback['radiation'], delta_age

    def prepare_new_state(self, old_state, math=np):
        new_state = old_state.copy()
        new_state.t += self.dt
        qc_sum = math.sum(old_state.qc)
        cooling_rate = bf.dynamic_cooling(self.w)
        dz = self.w * self.dt
        new_state.p, new_state.T, new_state.z = self.schemata['atmosphere'](old_state.p,
                                                          old_state.T,
                                                          np.mean(old_state.z),
                                                          dz)
        new_state.E = self.schemata['radiation'](new_state, self.microphysics)
        return new_state

def nucleation_slice(state, S_perturbations, microphysics):
    r_min = microphysics['r_min']
    particle_count = microphysics['particle_count']
    S = bf.relative_humidity(state.T, state.p, state.qv) - 1 + S_perturbations
    m = (S < bf.critical_super_saturation(r_min, state.T)) & (np.isclose(bf.radius(state.qc, particle_count, r_min), r_min))
    return np.invert(m)

def step(model, old_state, math=np):
    new_state = model.prepare_new_state(old_state, math=math)
    delta_Ts, delta_qvs, delta_qc, delta_T_E, delta_age = model.calculate_tendencies(new_state, math=math)
    new_state.qc = new_state.qc + delta_qc
    new_state.T += math.sum(delta_Ts) + delta_T_E
    new_state.qv += math.sum(delta_qvs)
    new_state.z = new_state.z + bf.fall_speed(bf.radius(new_state.qc, model.microphysics['particle_count'], model.microphysics['r_min'])) * model.dt
    new_state.age[delta_age<=0.] = 0.
    new_state.age = new_state.age + delta_age
    return new_state
