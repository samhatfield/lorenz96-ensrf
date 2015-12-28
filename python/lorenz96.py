import numpy as np
from math import sqrt

F = 8
DT = 0.05

def lorenz96(prev_state, model_error = False):
    # Number of X variables
    num_x = len(prev_state)

    k1 = florenz96(prev_state, num_x)
    k2 = florenz96(prev_state+0.5*DT*k1, num_x)
    k3 = florenz96(prev_state+0.5*DT*k2, num_x)
    k4 = florenz96(prev_state+DT*k3, num_x)

    next_state = prev_state + (DT/6)*(k1 + 2*k2 + 2*k3 + k4)

    # Include normally distributed errors, if specified
    if model_error:
        # Model error can be modified here
        sigma_x = np.ones(num_x)
        return next_state + sqrt(DT) * sigma_x * np.random.randn(num_x)
    else:
        return next_state

def florenz96(prev_state, K):
    return np.array([dXdT(prev_state, i, K) for i in range(K)])

def dXdT(X, k, K):
    return (X[(k+1)%K] - X[k-2])*X[k-1] - X[k] + F

