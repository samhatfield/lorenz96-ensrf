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

# Perform an assimilation every x timesteps
assim_freq = 1

#===============================================================================
# Spin up
#===============================================================================

print('Spin up...')

spins = 5000

# Initial perturbation
initial_truth = 8 * np.ones(num_x)
initial_truth[19] = 8.008

for _ in range(spins):
    initial_truth = lorenz96(initial_truth)

#===============================================================================
# Obtain truth run
#===============================================================================

print('Obtaining truth run...')

truth_run = [initial_truth]
for _ in range(num_steps-1):
    truth_run.append(lorenz96(truth_run[-1]))

#===============================================================================
# Extract and perturb observations
#===============================================================================

print('Obtaining observations...')

# For now, every X variable is observed
observations = np.array(truth_run)

var_o = 0.1 * np.ones(num_x)
sigma_o =  np.array([sqrt(x) for x in var_o])
Ro = np.diag(var_o)

# Perturb observations
observations = [ob + sigma_o * randn(num_x) for ob in observations]

#===============================================================================
# Setup filtering
#===============================================================================

# Ensemble size
num_ens = 50

# Initialization step
truth_mean = np.mean(truth_run, axis=0)
var0 = 3 * np.ones(num_x)
sigma0 = np.array([sqrt(x) for x in var0])

# Define ensemble and perturb members
print('Building initial ensemble...')

ensemble = [initial_truth + sigma0 * randn(num_x) for _ in range(num_ens)]

# Store first step
analysis_history = [np.mean(ensemble, axis=0)]

#===============================================================================
# 'Free' run; a single run initialised from the time-mean of the truth run, 
# for comparison
#===============================================================================

print('Obtaining free run...')

free_run = [truth_mean]
for _ in range(num_steps-1):
    free_run.append(lorenz96(free_run[-1]))

#===============================================================================
# Run filter
#===============================================================================

print('Running filter...')

for step in range(num_steps):
    if step % 10 == 0:
        print('Step %d' % step)

    # Forecast step
    ensemble = [lorenz96(member, model_error=True) for member in ensemble]

    # Time to perform an analysis?
    if step % assim_freq == 0:
        # Analysis step (for now, copy analysis from MATLAB version. Slow matrix
        # operations can be optimised out later)
        ensemble = assimilate(ensemble, observations[step], Ro)

    # Store ensemble mean
    analysis_history.append(np.mean(ensemble, axis=0))

#===============================================================================
# Plot results
#===============================================================================

# Plot comparison of truth with observations
plt.figure(1, facecolor='white')
truth_handle, = plt.plot([norm(x) for x in truth_run])
obs_handle, = plt.plot([norm(x) for x in observations])
plt.legend([truth_handle, obs_handle], ['Truth', 'Observations'])
plt.xlabel('MTUs')
plt.ylabel('Norm')

# Plot RMS differences between truth run and other runs
rms_free = []
rms_analy = []

for step in range(num_steps):
    normalization = 1.0/sqrt(num_x)

    rms_free.append(normalization * norm(truth_run[step] - free_run[step]))
    rms_analy.append(normalization * norm(truth_run[step] - analysis_history[step]))

fig = plt.figure(2, facecolor='white')
free_handle, = plt.plot(rms_free)
analy_handle, = plt.plot(rms_analy)
plt.legend([free_handle, analy_handle], ['Free', 'Analysis'])
plt.xlabel('MTUs')
plt.ylabel('Distance from truth state vector')

# Plot actual trajectories of forecast and truth norms
truth_norm = [norm(step) for step in truth_run]
analy_norm = [norm(step) for step in analysis_history]

fig = plt.figure(3, facecolor='white')
truth_handle, = plt.plot(truth_norm)
analy_handle, = plt.plot(analy_norm)
plt.legend([truth_handle, analy_handle], ['Truth', 'Analysis'])
plt.xlabel('MTUs')
plt.ylabel('State vector magnitude')

plt.show()
