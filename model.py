from scipy.integrate import odeint
import boxmodel_functions as bf
import boxmodel_constants as bc
from state import State

class Model(object):

    def __init__(self):
        self.r_min = (1e-6,)
        self.particle_count = (50e6,)
        self.rad_adjust = False
        self.T_env = 250.
        self.dt = 1.

    def run(self, logger):
        state = self.initial_state()
        logger.log_state(state)
        while not self.is_converged(state):
            state = self.step(state)
            logger.log_state(state)

    def initial_state(self):
        return State(0, 280, 100000, 0.010, (0,))

    def is_converged(self, state):
        return state.t >= 10 * 60

    def step(self, old_state):
        new_state = old_state.copy() 
        new_state.t += self.dt
        qc_sum = sum(old_state.qc)
        cooling_rate = bf.thermal_radiative_cooling_rate(old_state.T, qc_sum, self.T_env)
        new_state.T += cooling_rate * self.dt 

        delta_Ts, delta_qvs, new_state.qc = zip(*[self.condensation(old_state.T, 
                                                                    old_state.p, 
                                                                    old_state.qv, 
                                                                    qc_sum, 
                                                                    qc, 
                                                                    particle_count, 
                                                                    r_min)
                                                  for qc, particle_count, r_min
                                                  in zip(old_state.qc, self.particle_count, self.r_min)])
        new_state.T += sum(delta_Ts)
        new_state.qv += sum(delta_qvs)
        return new_state

    def condensation(self, T, p, qv, qc_sum, qc, particle_count, r_min):
        r_old = max(r_min, bf.radius(qc, particle_count))
      
        es = bf.saturation_pressure(T)
        S = bf.relative_humidity(T, p, qv) - 1
      
        if self.rad_adjust:
            E = bf.thermal_radiation(T, qc_sum)
        else:
            E = 0
       
        r_new = odeint(bf.differential_growth_by_condensation, r_old, [0, self.dt], args=(E, es, T, S), mxstep=2000)[1,0]
      
        qc_new = bf.cloud_water(particle_count, r_new)
        
        delta_qc = (qc_new - qc)
      
        delta_T = delta_qc * bc.H_LAT / bc.C_P
        delta_qv = -delta_qc
        return delta_T, delta_qv, qc_new
