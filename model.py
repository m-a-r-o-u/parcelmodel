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
        self.microphysics = {'r_min' : np.array(model_parameters['r_min']),
                             'particle_count' : np.array(model_parameters['particle_count'])}
        self.schemata = {'radiation': model_parameters['radiation_function'],
                         'turbulence': model_parameters['turbulence_schema'],
                         'atmosphere': model_parameters['atmosphere_schema']}
        self.w = model_parameters['w']
        self.dt = model_parameters['dt']
        self.t_max = model_parameters['t_max']
        initial_state = initial_state.copy()
        initial_state.qc = np.array(initial_state.qc, dtype='float')
        self._initial_state = initial_state
        self.output_step = model_parameters['output_step']
        assert len(self.microphysics['r_min']) == len(self.microphysics['particle_count'])
        assert len(self.microphysics['r_min']) == len(self._initial_state.qc)
        self.executer = executer(step, self, initial_state, int(self.output_step / float(self.dt)))
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

        delta_Ts_act, delta_qvs_act, delta_qc_act = np.zeros((3, len(state.qc)))
#        S = bf.relative_humidity(state.T, state.p, state.qv) - 1 + S_perturbations
#        m_act = np.isclose(state.qc, 0) & (S > bf.critical_super_saturation(self.microphysics['r_min'], state.T))
#
#        r_act = 2. / 3. * bf.kelvins_parameter(state.T) / bf.critical_super_saturation(self.microphysics['r_min'][m_act], state.T)
#        delta_qc_act[m_act] = bf.cloud_water(self.microphysics['particle_count'][m_act], 
#                r_act,
#                self.microphysics['r_min'][m_act])
#        delta_qvs_act[m_act] = -delta_qc_act[m_act]
#        import boxmodel_constants as c
#        delta_Ts_act[m_act] = delta_qc_act[m_act] * c.H_LAT / c.C_P
#
#        #if np.any(m_act) : print m_act, np.max(r_act), np.mean(self.microphysics['r_min']), state.t


        def condensation(qc, particle_count, r_min, S_perturbation, E):
            return bf.condensation(state.T, state.p, state.qv, qc_sum, qc, particle_count, r_min, self.dt, E, S_perturbation, math=math)
        m = nucleation_slice(state, S_perturbations, self.microphysics)
        delta_Ts, delta_qvs, delta_qc = np.zeros((3, len(state.qc)))
        delta_Ts[m], delta_qvs[m], delta_qc[m] = condensation(state.qc[m],
                                                              self.microphysics['particle_count'][m],
                                                              self.microphysics['r_min'][m],
                                                              S_perturbations[m],
                                                              state.E[0])
        return delta_Ts + delta_Ts_act, delta_qvs + delta_qvs_act, delta_qc + delta_qc_act

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
    delta_Ts, delta_qvs, delta_qc = model.calculate_tendencies(new_state, math=math)
    new_state.qc = new_state.qc + delta_qc
    new_state.T += math.sum(delta_Ts)
    new_state.qv += math.sum(delta_qvs)
    new_state.z = new_state.z + bf.fall_speed(bf.radius(new_state.qc, model.microphysics['particle_count'], model.microphysics['r_min'])) * model.dt
    return new_state
