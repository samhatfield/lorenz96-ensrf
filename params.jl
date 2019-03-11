# Frequency with which to write output
const write_freq = 1

# Output filename
const outfile = "output.nc"

# Time step
const Δt = 0.005

# Start and end times, and number of time steps
const fin = 2.0
const nsteps = Int(fin / Δt)

# State vector properties
const nx = 8
const ny = 8
const nz = 8

# State dimension of full model
const truthdim = nx + nx*ny + nx*ny*nz

# State dimension of parametrised model
const statedim = nx + nx*ny

# Number of ensemble members
const nens = 400

# Frequency of assimilations, i.e. 1 = every timestep
const assim_freq = 10

# Observation error variance
const σ² = 0.1

# Dimension of observation vector
const obsdim = nx*ny

# Covariance inflation factor
const ρ = 1.15