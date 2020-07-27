import numpy as np
from .state import State


class NormalExecuter(object):
    def __init__(self, function, model, state, n_steps):
        self.func = function
        self.model = model
        self.n_steps = n_steps

    def step(self, state):
        for _ in range(self.n_steps):
            state = self.func(self.model, state, math=np)
        return state


class TheanoExecuter():
    def __init__(self, function, model, state, n_steps):
        import theano.tensor as TT
        import theano
        variables = [TT.dscalar('t'),
                     TT.dscalar('T'),
                     TT.dscalar('p'),
                     TT.dscalar('qv'),
                     TT.dvector('qc')]

        def step_wrapper(t, T, p, qv, qc):
            state = State(t, T, p, qv, qc)
            new_state = function(model, state, math=TT)
            return new_state.t, new_state.T, new_state.p, new_state.qv, new_state.qc

        result, updates = theano.scan(fn=step_wrapper,
                                      outputs_info=variables,
                                      n_steps=n_steps)

        final_result = [res[-1] for res in result]

        self.func = theano.function(
            inputs=variables, outputs=final_result, updates=updates)

    def step(self, state):
        t, T, p, qv, qc = self.func(state.t,
                                    state.T,
                                    state.p,
                                    state.qv,
                                    state.qc)
        return State(t, T, p, qv, qc)


DISPATCHER = {'numpy': NormalExecuter,
              'theano': TheanoExecuter}
