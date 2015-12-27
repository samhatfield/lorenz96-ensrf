from numpy.random import randn
from numpy.linalg import norm
import numpy as np
from math import sqrt
from lorenz96 import lorenz96, DT
import matplotlib.pyplot as plt

#===============================================================================
# Setup
#===============================================================================

# Number of X variables
num_x = 40

# Starting and ending times
t1, t2 = DT, 6

# Set total number of model steps (last step isn't included by 
# arange, so it is added separately
num_steps = len(np.arange(t1, t2, DT)) + 1

#===============================================================================
# Spin up
#===============================================================================

spins = 5000
state = 8 * np.ones(num_x)
state[19] = 8.008

for _ in range(spins):
    state = lorenz96(state)

initial_truth = state 

#===============================================================================
# Obtain truth run
#===============================================================================

truth = []
next_x = initial_truth
for _ in range(num_steps):
    truth.append(next_x)
    next_x = lorenz96(next_x)

#===============================================================================
# Extract and perturb observations
#===============================================================================

# For now, every X variable is observed
obs = np.array(truth)

# Set frequency of observations (same as DT for now) (must be multiple of DT)
DT_obs = DT
assert(DT_obs % DT == 0.0)

var_o = 0.01 * np.ones(num_x)
sigma_o = np.array([sqrt(x) for x in var_o])
# r_o = 

# Perturb observations
for i in range(len(obs)):
    obs[i] = obs[i] + sigma_o * randn(num_x)

#===============================================================================
# Setup filtering
#===============================================================================

# Ensemble size
num_ens = 50

# Number of forecast steps between each analysis step
num_fore_steps = int(DT_obs/DT)

# Total number of filter steps
num_filt_steps = len(np.arange(t1, t2, DT_obs))

# Set inflation factor
rho = 1.4

# Initialization step
Xa0 = np.mean(truth, axis=0)
var0 = 3 * np.ones(num_x)
sigma0 = np.array([sqrt(x) for x in var0])

# Define ensemble and perturb members
ensemble = []
for _ in range(num_ens):
    ensemble.append(Xa0 + sigma0 * randn(num_x))

# Store first step
forecast_history = analysis_history = [np.mean(ensemble, axis=0)]

#===============================================================================
# Run filter
#===============================================================================

for filt_step in range(num_filt_steps):
    print 'Step %d' % filt_step

    # Forecast steps
    for _ in range(num_fore_steps):
        ensemble = [lorenz96(member) for member in ensemble]


    ens_mean = np.mean(ensemble, axis=0)

    # Save ensemble state prior to analysis
    forecast_history.append(ens_mean)

    print 'End forecast'

#===============================================================================
# Plot results
#===============================================================================

rms_fore = []
for filt_step in range(num_filt_steps):
    skip = filt_step * num_fore_steps

    rms_fore.append(1.0/sqrt(num_x) * norm(truth[skip] - forecast_history[skip]))

# Plot RMS differences between truth and other runs
plt.plot(rms_fore)
plt.show()
