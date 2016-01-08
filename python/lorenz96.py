import numpy as np

F = 20
DT = 0.05
H = 1
C = 4
B = 10

"""
Runge-Kutta integration scheme.
"""
def lorenz96(prev_state, size):
    X = prev_state[:size[0]]
    Y = prev_state[size[0]:]

    k1 = dXdT(X, Y)
    l1 = dYdT(X, Y)
    k2 = dXdT(X+0.5*DT*k1, Y+0.5*DT*l1)
    l2 = dYdT(X+0.5*DT*k1, Y+0.5*DT*l1)
    k3 = dXdT(X+0.5*DT*k2, Y+0.5*DT*l2)
    l3 = dYdT(X+0.5*DT*k2, Y+0.5*DT*l2)
    k4 = dXdT(X+DT*k3, Y+DT*l3)
    l4 = dYdT(X+DT*k3, Y+DT*l3)

    X = X + (DT/6)*(k1 + 2*k2 + 2*k3 + k4)
    Y = Y + (DT/6)*(l1 + 2*l2 + 2*l3 + l4)

    return np.concatenate([X, Y])

def dXdT(X, Y):
    sum_y = np.sum(np.split(Y, X.shape[0]), axis=1)
    return (np.roll(X, -1)-np.roll(X, 2))*np.roll(X, 1) - X + F - (H*C/B)*sum_y

def dYdT(X, Y):
    X_rpt = np.repeat(X, Y.shape[0]/X.shape[0])
    return (np.roll(Y, 1)-np.roll(Y, -2))*C*B*np.roll(Y, -1)-C*Y+(H*C/B)*X_rpt
