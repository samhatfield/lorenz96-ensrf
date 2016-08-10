# Start and end times, and number of time steps
dt = 0.005
fin = 12.0
n_steps = int(fin / dt)

n_x = 8
n_y = 32

# State vector dimension
state_dim = n_x + n_x*n_y

# Number of ensemble members
n_ens = 500

assim_freq = 1

# Observation error variance
y_var = 0.1

# Only observe Y variables
obs_dim = n_x*n_y
