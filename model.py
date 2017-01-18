from scipy.integrate import odeint
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
        self.r_min = model_parameters['r_min']
        self.particle_count = model_parameters['particle_count']
        self.radiation = model_parameters.get('radiation', False)
        self.T_env = model_parameters['T']
        self.dt = model_parameters['dt']
        self.t_max = model_parameters['t_max']
        self._initial_state = initial_state
        assert len(self.r_min) == len(self.particle_count)
        assert len(self.r_min) == len(self._initial_state.qc)

    def run(self, logger):
        state = self.initial_state()
        logger.set_units(self.units)
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
        delta_Ts, delta_qvs, new_state.qc = self.calculate_tendencies(new_state)

        new_state.T += sum(delta_Ts)
        new_state.qv += sum(delta_qvs)
        return new_state

    def calculate_tendencies(self, state):
        qc_sum = sum(state.qc)
        def condensation(qc, particle_count, r_min):
            return self.condensation(state.T, state.p, state.qv, qc_sum, qc, particle_count, r_min)
        delta_Ts, delta_qvs, state.qc = zip(*map(condensation, state.qc, self.particle_count, self.r_min))
        return delta_Ts, delta_qvs, state.qc

    def prepare_new_state(self, old_state):
        new_state = old_state.copy()
        new_state.t += self.dt
        qc_sum = sum(old_state.qc)
        cooling_rate = bf.thermal_radiative_cooling_rate(old_state.T, qc_sum, self.T_env)
        new_state.T += cooling_rate * self.dt
        return new_state

    def condensation(self, T, p, qv, qc_sum, qc, particle_count, r_min):
        r_old = max(r_min, bf.radius(qc, particle_count))
      
        es = bf.saturation_pressure(T)
        S = bf.relative_humidity(T, p, qv) - 1
      
        if self.radiation:
            E = bf.thermal_radiation(T, qc_sum)
        else:
            E = 0
       
        r_new = odeint(bf.differential_growth_by_condensation, r_old, [0, self.dt], args=(E, es, T, S), mxstep=2000)[1,0]
      
        qc_new = bf.cloud_water(particle_count, r_new)
        
        delta_qc = (qc_new - qc)
      
        delta_T = delta_qc * bc.H_LAT / bc.C_P
        delta_qv = -delta_qc
        return delta_T, delta_qv, qc_new
