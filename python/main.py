"""
Authored by Sam Hatfield, heavily based on MATLAB code by Aneesh
Subramanian.
"""

from numpy.random import randn
from numpy.linalg import norm
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt

from lorenz96 import lorenz96, DT
from analysis import assimilate

#===============================================================================
# Setup
#===============================================================================

# Number of X variables
num_x = 40

# Starting and ending times
t1, t2 = DT, 12 

# Set total number of model steps (last step isn't included by 
# arange, so it is added separately
num_steps = len(np.arange(t1, t2, DT)) + 1

#===============================================================================
# Spin up
#===============================================================================

spins = 5000

# Initial perturbation (Todo: why is this necessary?)
initial_truth = 8 * np.ones(num_x)
initial_truth[19] = 8.008

for _ in range(spins):
    initial_truth = lorenz96(initial_truth)

#===============================================================================
# Obtain truth run
#===============================================================================

truth_run = [initial_truth]
for _ in range(num_steps-1):
    truth_run.append(lorenz96(truth_run[-1]))

#===============================================================================
# Extract and perturb observations
#===============================================================================

# For now, every X variable is observed
obs = np.array(truth_run)

var_o = 0.1 * np.ones(num_x)
sigma_o =  np.array([sqrt(x) for x in var_o])
Ro = np.diag(var_o)

# Perturb observations
obs = [ob + sigma_o * randn(num_x) for ob in obs]

#===============================================================================
# Setup filtering
#===============================================================================

# Ensemble size
num_ens = 50

# Initialization step
Xa0 = np.mean(truth_run, axis=0)
var0 = 3 * np.ones(num_x)
sigma0 = np.array([sqrt(x) for x in var0])

# Define ensemble and perturb members
ensemble = [Xa0 + sigma0 * randn(num_x) for _ in range(num_ens)]

# Store first step
forecast_history = [np.mean(ensemble, axis=0)]
analysis_history = [np.mean(ensemble, axis=0)]

#===============================================================================
# 'Free' run; a single run initialised from the time-mean of the truth run, 
# for comparison
#===============================================================================

free_run = [Xa0]
for _ in range(num_steps-1):
    free_run.append(lorenz96(free_run[-1]))

#===============================================================================
# Run filter
#===============================================================================

for step in range(num_steps-1):
    print 'Step %d' % step

    # Forecast step
    ensemble = [lorenz96(member) for member in ensemble]

    ens_mean = np.mean(ensemble, axis=0)

    # Save ensemble state prior to analysis
    forecast_history.append(ens_mean)

    print 'End forecast'

    # Analysis step (for now, copy analysis from MATLAB version. Slow matrix
    # operations can be optimised out later)
    ensemble = assimilate(ensemble, obs[step], Ro)

    analysis_history.append(np.mean(ensemble, axis=0))

    print 'End analysis'

#===============================================================================
# Plot results
#===============================================================================

rms_free = []
rms_fore = []
rms_anal = []

for step in range(num_steps):
    normalization = 1.0/sqrt(num_x)

    rms_free.append(normalization * norm(truth_run[step] - free_run[step]))
    rms_fore.append(normalization * norm(truth_run[step] - forecast_history[step]))
    rms_anal.append(normalization * norm(truth_run[step] - analysis_history[step]))

# Plot RMS differences between truth run and other runs
free_handle, = plt.plot(rms_free)
fore_handle, = plt.plot(rms_fore)
anal_handle, = plt.plot(rms_anal)
plt.legend([free_handle, fore_handle, anal_handle], ['Free', 'Forecast', 'Analysis'])
plt.show()

