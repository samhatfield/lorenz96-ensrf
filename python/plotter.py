import json
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import sqrt
import numpy as np
from sys import argv
from matplotlib import rc
import re
import matplotlib as mpl

mpl.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})
mpl.rc('text', usetex=True)
plt.style.use('ggplot')
mpl.rc('text.latex', preamble='\\usepackage{sfmath}')

with open(argv[1]) as data_file:
    results = json.load(data_file)

filter_history = results['filter_history']
truth_run = results['truth_run']
num_steps = results['params']['num_steps']
N_X = results['params']['N_X']
N_Y = results['params']['N_Y']
assim_freq = results['params']['assim_freq']

# Plot actual trajectories of filter and truth norms, along with ensemble
# spread
truth_norm = [norm(step[:N_X]) for step in truth_run]
filt_norm = np.array([norm(step['x_mean_norm']) for step in filter_history[1:]])
x_std_norm = np.array([step['x_std_norm'] for step in filter_history[1:]])

fig = plt.figure(figsize=(12,6))
truth_handle, = plt.plot(truth_norm)
filt_handle, = plt.plot(filt_norm)
plt.fill_between(range(len(filt_norm)), filt_norm-x_std_norm, filt_norm+x_std_norm, facecolor='green', alpha=0.5, edgecolor='none')
leg = plt.legend([truth_handle, filt_handle ], ['Truth', 'Analysis'], frameon=True)
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
