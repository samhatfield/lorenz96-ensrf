import numpy as np
from math import sqrt

F = 8
DT = 0.05

"""
Runge-Kutta integration scheme. This operates on the entire state vector at the
same time.
"""
def lorenz96(prev_state, model_error = False):
    num_x = len(prev_state)

    k1 = dXdT(prev_state)
    k2 = dXdT(prev_state+0.5*DT*k1)
    k3 = dXdT(prev_state+0.5*DT*k2)
    k4 = dXdT(prev_state+DT*k3)

    next_state = prev_state + (DT/6)*(k1 + 2*k2 + 2*k3 + k4)

    # Include normally distributed errors, if specified
    if model_error:
        # Model error can be modified here
        return next_state + sqrt(DT) * np.random.randn(num_x)
    else:
        return next_state

def dXdT(X):
    return (np.roll(X, -1) - np.roll(X, 2)) * np.roll(X, 1) - X + F

