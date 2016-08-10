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
from analysis import ensrf_assimilate
from observation import observe
from params import *

from datetime import datetime

#===============================================================================
# Setup
#===============================================================================

# Keep track of execution time
start = datetime.now()

# Set the seed
seed(0)

#===============================================================================
# Spin up
#===============================================================================

print('Spin up...')

spins = 5000

# Build state vector and perturb it
initial_truth = np.concatenate([8 * np.ones(n_x), 0.5*np.ones(n_x*n_y)])
initial_truth[3] = 8.008

for _ in range(spins):
    initial_truth = lorenz96(initial_truth, (n_x, n_y))

#===============================================================================
# Obtain truth run
#===============================================================================

print('Generating truth...')

truth_run = [initial_truth]
for _ in range(n_steps-1):
    truth_run.append(lorenz96(truth_run[-1], (n_x, n_y)))

#===============================================================================
# Extract and perturb observations
#===============================================================================

print('Extracting observations...')

# Make observations
observations = observe(np.transpose(truth_run))

# Define observational error covariance matrix (diagonal matrix of variances)
obs_covar = y_var * np.identity(obs_dim)

# Perturb observations
observations = [ob + sqrt(y_var) * randn(num_obs) for ob in np.transpose(observations)]

#===============================================================================
# Setup filtering
#===============================================================================

print('Generating ensemble...')

# Get climatology from time average of truth
climatology_mean = np.mean(truth_run, axis=0)
climatology_std = np.std(truth_run, axis=0)

ensemble = np.empty((n_x + n_x*n_y, num_ens))
for i in range(num_ens):
    ensemble[:,i] = climatology_mean + climatology_std * randn(state_dim)

#===============================================================================
# Run filter
#===============================================================================

print('Running filter...')

# We skip the first step because the initial conditions are the first step
for step in range(1, n_steps):
    if step % 10 == 0:
        print('Step %d' % step)

    # Forecast step
    for i in range(num_ens):
        ensemble[:,i] = lorenz96(ensemble[:,i], (n_x, n_y))

    # Time to perform an analysis?
    if step % assim_freq == 0:
        # Analysis step
        ensemble = ensrf_assimilate(ensemble, observations[step], Ro, sigma_o)

    # Store ensemble mean
    filter_history.append({
        'mean': np.mean(ensemble, axis=1),
        'x_mean_norm': np.mean(norm(ensemble[:n_x,:], axis=0)),
        'x_std_norm': np.std(norm(ensemble[:n_x,:], axis=0)),
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

# Collect results
results = {
    'time_elapsed': str(datetime.now() - start),
    'params': {
        'n_x': n_x,
        'n_y': n_y,
        'DT': DT,
        'num_ens': num_ens,
        'assim_freq': assim_freq,
        'n_steps': n_steps
    },
    'filter_history': filter_history,
    'truth_run': truth_run,
    'observations': observations
}

with open('results.json', 'w') as json_file:
    json.dump(results, json_file)
