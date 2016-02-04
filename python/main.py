"""
Authored by Sam Hatfield, heavily based on MATLAB code by Aneesh
Subramanian.
"""

from numpy.random import randn, seed
from numpy.linalg import norm
import numpy as np
from math import sqrt
import json

from lorenz96 import lorenz96, DT, H, C, B
from analysis import assimilate
from observation import observe
from params import N_X, N_Y

from datetime import datetime

#===============================================================================
# Setup
#===============================================================================

# Keep track of execution time
start = datetime.now()

# Starting and ending times
t1, t2 = 0, 25

# Set total number of model steps (last step isn't included by 
# arange, so it is added separately
num_steps = len(np.arange(t1, t2, DT)) + 1

# Perform an assimilation every x timesteps
assim_freq = 1

# Set the seed
seed(0)

#===============================================================================
# Spin up
#===============================================================================

print('Spin up...')

spins = 5000

# Build state vector and perturb it
initial_truth = np.concatenate([8 * np.ones(N_X), 0.5*np.ones(N_X*N_Y)])
initial_truth[3] = 8.008

for _ in range(spins):
    initial_truth = lorenz96(initial_truth, (N_X, N_Y))

#===============================================================================
# Obtain truth run
#===============================================================================

print('Obtaining truth run...')

truth_run = [initial_truth]
for _ in range(num_steps-1):
    truth_run.append(lorenz96(truth_run[-1], (N_X, N_Y)))

#===============================================================================
# Extract and perturb observations
#===============================================================================

print('Obtaining observations...')

# Make observations
observations = observe(np.transpose(truth_run))

# Dimension of observed state vector
num_obs = observations.shape[0]

var_o = 0.1 * np.ones(num_obs)
sigma_o =  np.array([sqrt(x) for x in var_o])
Ro = np.diag(var_o)

# Perturb observations
observations = [ob + sigma_o * randn(num_obs) for ob in np.transpose(observations)]

#===============================================================================
# Setup filtering
#===============================================================================

# Ensemble size
num_ens = 500

# Initialization step
var0 = 3 * np.ones(N_X + N_X*N_Y)
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
ensemble = np.empty((N_X + N_X*N_Y, num_ens))
for i in range(num_ens):
    ensemble[:,i] = initial_truth + sigma0 * randn(N_X + N_X*N_Y)

# Store first step
filter_history = [{
    'mean': np.mean(ensemble, axis=1),
    'x_mean_norm': np.mean(norm(ensemble[:N_X,:], axis=0)),
    'x_upp_norm': np.max(norm(ensemble[:N_X,:], axis=0)),
    'x_low_norm': np.min(norm(ensemble[:N_X,:], axis=0))
}]

#===============================================================================
# 'Free' run; a single run initialised from the first ensemble member, 
# without analysis for comparison
#===============================================================================

print('Obtaining free run...')

free_run = [ensemble[:,0]]
for _ in range(num_steps-1):
    free_run.append(lorenz96(free_run[-1], (N_X, N_Y)))

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
        ensemble[:,i] = lorenz96(ensemble[:,i], (N_X, N_Y))

    # Time to perform an analysis?
    if step % assim_freq == 0:
        # Analysis step
        ensemble = assimilate(ensemble, observations[step], Ro, sigma_o)

    # Store ensemble mean
    filter_history.append({
        'mean': np.mean(ensemble, axis=1),
        'x_mean_norm': np.mean(norm(ensemble[:N_X,:], axis=0)),
        'x_upp_norm': np.max(norm(ensemble[:N_X,:], axis=0)),
        'x_low_norm': np.min(norm(ensemble[:N_X,:], axis=0))
    })

print('Time elapsed: ' + str(datetime.now() - start))

#===============================================================================
# Save to JSON
#===============================================================================

# Make filter history JSON-parseable
for step in filter_history:
    for key, value in step.items():
        if key == 'mean':
            step['mean'] = step['mean'].tolist()

# Make truth run JSON-parseable
truth_run = [step.tolist() for step in truth_run]

# Make observations JSON-parseable
observations = [step.tolist() for step in observations]

# Make free run JSON-parseable
free_run = [step.tolist() for step in free_run]

# Collect results
results = {
    'time_elapsed': str(datetime.now() - start),
    'params': {
        'N_X': N_X,
        'N_Y': N_Y,
        'DT': DT,
        'num_ens': num_ens,
        'assim_freq': assim_freq,
        'num_steps': num_steps
    },
    'filter_history': filter_history,
    'truth_run': truth_run,
    'observations': observations,
    'free_run': free_run
}

with open('results.json', 'w') as json_file:
    json.dump(results, json_file)
