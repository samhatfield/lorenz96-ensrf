import matplotlib.pyplot as plt
import matplotlib as mpl
from sys import argv
from matplotlib import rc
import re
import yaml
from math import log
import numpy as np

mpl.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Arial']})
mpl.rc('text', usetex=True)
plt.style.use('seaborn-colorblind')
mpl.rc('text.latex', preamble='\\usepackage{sfmath}')

# Open data file
with open(argv[1], 'r') as stream:
    docs = stream.read().split('---')

    # Get file metadata
    meta = yaml.load(docs[0])

    # Get actual data
    data = docs[1][1:]

x_avg_norm = []
x_norm_std = []
truth = []
for row in data.splitlines():
    x_avg_norm.append(float(row.split()[0]))
    x_norm_std.append(float(row.split()[1]))
    truth.append(float(row.split()[2]))

# Process metadata
num_steps = int(meta['fin']/meta['dt'])
N_X = meta['n_x']
N_Y = meta['n_y']
N_Z = meta['n_z']
assim_freq = meta['assim_freq']
n_ens = meta['n_ens']

# Plot actual trajectories of filter and truth norms, along with ensemble
# spread
fig = plt.figure(3, facecolor='white', figsize=(12,6))
truth_handle, = plt.plot(truth)
filt_handle, = plt.plot(x_avg_norm)
x_avg_norm = np.array(x_avg_norm)
x_norm_std = np.array(x_norm_std)
plt.fill_between(range(len(x_avg_norm)), x_avg_norm-x_norm_std, x_avg_norm+x_norm_std, facecolor='green', alpha=0.5, edgecolor='none')
leg = plt.legend([truth_handle, filt_handle], ['Truth', 'Filter'], frameon=True)
rect = leg.get_frame()
rect.set_linewidth(0.0)
rect.set_alpha(0.7)

plt.xlabel('Time steps')
plt.ylabel('X norm')
title_timesteps = 'time step' if assim_freq == 1 else str(assim_freq) + ' time steps'
title_timesteps += ' (%d, %d, %d)' % (N_X, N_Y, N_Z)
title_timesteps += ' %d members' % n_ens
#plt.title('Assimilations every ' + title_timesteps)

# Get figure ID from file name
#p = re.compile(r'\S+/(\S+).yml')
#match = p.match(argv[1])
#figure_name = match.group(1)
#
## Save to PDF
#plt.savefig('figures/trajectories/%s.pdf' % figure_name)
#
#print('Output saved to figures/trajectories/%s.pdf' % figure_name)

plt.show()
