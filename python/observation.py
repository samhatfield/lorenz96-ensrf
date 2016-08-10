from params import n_x

"""
Defines our observation operation.
"""
def observe(state):
    if len(state.shape) == 1:
        return state[N_X:]
    else:
        return state[N_X:,:]
