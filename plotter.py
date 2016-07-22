import matplotlib.pyplot as plt
from matplotlib import rc
from iris import load, analysis, FUTURE
import iris.quickplot as qplt
import iris.plot as iplt

# Silence future error message
FUTURE.netcdf_promote = True

# Set style
rc('text', usetex=True)
plt.style.use('ggplot')
rc('text.latex', preamble='\\usepackage{sfmath}')

# Load cubes
ens, truth = tuple(load('output.nc', ['ensx', 'truthx']))

# Compute ensemble mean of x mean of ensemble, and x mean of truth
ens_x_mean = ens.collapsed(['x'], analysis.MEAN)
ens_mean = ens_x_mean.collapsed(['ens'], analysis.MEAN)
ens_std = ens_x_mean.collapsed(['ens'], analysis.STD_DEV)
truth_mean = truth.collapsed(['x'], analysis.MEAN)

# Get time coordinate
time = ens.coord('time')

# Plot truth, ens mean and ens spread
fig = plt.figure(figsize=(12,6))
ens_h, = iplt.plot(ens_mean)
truth_h, = iplt.plot(truth_mean)
plt.fill_between(
        time.points, (ens_mean-ens_std).data, (ens_mean+ens_std).data,
        facecolor=ens_h.get_color(), alpha=0.5, edgecolor='none'
)

# Add legend
leg = plt.legend([truth_h, ens_h], ['Truth', 'Ensemble'], frameon=True)
rect = leg.get_frame()
rect.set_linewidth(0.0)
rect.set_alpha(0.7)

plt.xlabel('MTUs')
plt.ylabel('$\\bar{X}$')

# Save to PDF
plt.savefig('output.pdf')

plt.show()
