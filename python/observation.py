from params import N_X

"""
Defines our observation operation.
"""
def observe(state):
    return state[N_X:,:]
