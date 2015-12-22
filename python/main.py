from numpy import arange, zeros, ones, array
from numpy.random import randn
from math import sqrt
from lorenz96 import lorenz96, DT

# Number of X variables
num_x = 40

# Starting and ending times
t1, t2 = DT, 6

# Set total number of model steps (last step isn't included by 
# arange, so it is added separately
num_steps = len(arange(t1, t2, DT)) + 1

# Initialize solution arrays
#truth = free = filter_a = filter_f = zeros((num_steps, num_x))
truth = []

# Spin up
spins = 5000
state = 8 * ones(num_x)
state[19] = 8.008

for _ in range(spins):
    state = lorenz96(state)

initial_truth = state 

# Set 'model error'
var_x = ones(num_x)
sigma_x = array([sqrt(x) for x in var_x])
alpha = 0.0

# Obtain truth run
next_x = initial_truth
for _ in range(num_steps):
    truth.append(next_x)
    next_x = lorenz96(next_x) + alpha*sqrt(DT)*sigma_x*randn(num_x)

print truth[-1]
