from numpy import array

F = 8
DT = 0.05

def lorenz96(prev_state):
    k1 = florenz96(prev_state)
    k2 = florenz96(prev_state+0.5*DT*k1)
    k3 = florenz96(prev_state+0.5*DT*k2)
    k4 = florenz96(prev_state+DT*k3)

    return prev_state + (DT/6)*(k1 + 2*k2 + 2*k3 + k4)

def florenz96(prev_state):
    num_x = len(prev_state)

    return array([dXdT(prev_state, i) for i in range(num_x)])

def dXdT(X, k):
    K = len(X)
    return (X[(k+1)%K] - X[k-2])*X[k-1] - X[k] + F

