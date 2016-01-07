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
t1, t2 = DT, 24 

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
num_ens = 500

# Initialization step
var0 = 3 * np.ones(num_x)
sigma0 = np.array([sqrt(x) for x in var0])

# Define ensemble and perturb members
print('Building initial ensemble...')

# Ensemble members are arranged like this...
# [F S T     L
#  i e h     a
#  r c i ... s
#  s o r     t
#  t n d
#    d        ]
ensemble = np.empty((num_x, num_ens))
for i in range(num_ens):
    ensemble[:,i] = initial_truth + sigma0 * randn(num_x)

# Store first step
filter_history = [{
    'mean': np.mean(ensemble, axis=1),
    'mean_norm': np.mean(norm(ensemble, axis=0)),
    'upp_norm': np.max(norm(ensemble, axis=0)),
    'low_norm': np.min(norm(ensemble, axis=0))
}]

#===============================================================================
# 'Free' run; a single run initialised from the time-mean of the truth run, 
# for comparison
#===============================================================================

print('Obtaining free run...')

free_run = [ensemble[:,0]]
for _ in range(num_steps-1):
    free_run.append(lorenz96(free_run[-1]))

#===============================================================================
# Run filter
#===============================================================================

print('Running filter...')

# We skip the first step because the initial conditions are the first step
for step in range(1, num_steps):
    if step % 10 == 0:
        print('Step %d' % step)

    # Forecast step
    for i in range(num_ens):
        ensemble[:,i] = lorenz96(ensemble[:,i])

    # Time to perform an analysis?
    if step % assim_freq == 0:
        # Analysis step
        ensemble = assimilate(ensemble, observations[step], Ro, sigma_o)

    # Store ensemble mean
    filter_history.append({
        'mean': np.mean(ensemble, axis=1),
        'mean_norm': np.mean(norm(ensemble, axis=0)),
        'upp_norm': np.max(norm(ensemble, axis=0)),
        'low_norm': np.min(norm(ensemble, axis=0))
    })

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
rms_filt = []

for step in range(num_steps):
    normalization = 1.0/sqrt(num_x)

    rms_free.append(normalization * norm(truth_run[step] - free_run[step]))
    rms_filt.append(normalization * norm(truth_run[step] - filter_history[step]['mean']))

fig = plt.figure(2, facecolor='white')
free_handle, = plt.plot(rms_free)
filt_handle, = plt.plot(rms_filt)
plt.legend([free_handle, filt_handle], ['Free', 'Analysis'])
plt.xlabel('MTUs')
plt.ylabel('Distance from truth state vector')

# Plot actual trajectories of filter, free and truth norms, along with ensemble
# spread
truth_norm = [norm(step) for step in truth_run]
filt_norm = [norm(step['mean_norm']) for step in filter_history]
free_norm = [norm(step) for step in free_run]
upp_norm = [step['upp_norm'] for step in filter_history]
low_norm = [step['low_norm'] for step in filter_history]

fig = plt.figure(3, facecolor='white')
truth_handle, = plt.plot(truth_norm)
filt_handle, = plt.plot(filt_norm)
free_handle,  = plt.plot(free_norm, color=(1,0.6,0.6,0.5))
plt.fill_between(range(num_steps), low_norm, upp_norm, facecolor='green', alpha=0.5, edgecolor='none')
plt.legend([truth_handle, filt_handle, free_handle], ['Truth', 'Analysis', 'Free'])
plt.xlabel('MTUs')
plt.ylabel('State vector magnitude')

plt.show()
