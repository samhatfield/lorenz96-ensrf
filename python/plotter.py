import json
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import sqrt
import numpy as np
from sys import argv
from matplotlib import rc
import re

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

with open(argv[1]) as data_file:
    results = json.load(data_file)

filter_history = results['filter_history']
truth_run = results['truth_run']
observations = results['observations']
free_run = results['free_run']
num_steps = results['params']['num_steps']
N_X = results['params']['N_X']
N_Y = results['params']['N_Y']
assim_freq = results['params']['assim_freq']

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
    normalization = 1.0/sqrt(N_X)

    rms_free.append(normalization * norm(np.array(truth_run[step]) - np.array(free_run[step])))
    rms_filt.append(normalization * norm(np.array(truth_run[step]) - np.array(filter_history[step]['mean'])))

fig = plt.figure(2, facecolor='white')
free_handle, = plt.plot(rms_free)
filt_handle, = plt.plot(rms_filt)
plt.legend([free_handle, filt_handle], ['Free', 'Analysis'])
plt.xlabel('MTUs')
plt.ylabel('Distance from truth state vector')

# Plot actual trajectories of filter, free and truth norms, along with ensemble
# spread
truth_norm = [norm(step[:N_X]) for step in truth_run]
filt_norm = [norm(step['x_mean_norm']) for step in filter_history]
free_norm = [norm(step[:N_X]) for step in free_run]
x_upp_norm = [step['x_upp_norm'] for step in filter_history]
x_low_norm = [step['x_low_norm'] for step in filter_history]

fig = plt.figure(3, facecolor='white', figsize=(12,6))
truth_handle, = plt.plot(truth_norm)
filt_handle, = plt.plot(filt_norm)
free_handle,  = plt.plot(free_norm, color=(1,0.6,0.6,0.5))
plt.fill_between(range(num_steps), x_low_norm, x_upp_norm, facecolor='green', alpha=0.5, edgecolor='none')
leg = plt.legend([truth_handle, filt_handle, free_handle], ['Truth', 'Analysis', 'Free'], frameon=True)
rect = leg.get_frame()
rect.set_linewidth(0.0)
rect.set_alpha(0.7)

plt.xlabel('MTUs')
plt.ylabel('State vector magnitude')
title_timesteps = 'MTU' if assim_freq == 1 else str(assim_freq) + ' MTUs'
title_timesteps += ' (%d, %d)' % (N_X, N_Y)
plt.title('Assimilations every ' + title_timesteps)
plt.ylim([10, 35])


p = re.compile(r'\S+/(\S+).json')
match = p.match(argv[1])
figure_name = match.group(1)

plt.savefig('figures/%s.pdf' % figure_name)

plt.show()
