#-------------------------------------------------------------------------------
# %% Import all function definitions and parameters
#-------------------------------------------------------------------------------

include("setup.jl")
include("params.jl")
include("lorenz96.jl")
include("observation.jl")
include("io.jl")
include("analysis.jl")

using LinearAlgebra: I
using Random: randn

#-------------------------------------------------------------------------------
# %% Spin up
#-------------------------------------------------------------------------------

println("Spinning up...")

truth_full = spinup()

#-------------------------------------------------------------------------------
# %% Truth run
#-------------------------------------------------------------------------------

println("Generating truth...")

truth = zeros(statedim, nsteps)
truth[:,1] = truth_full[1:nx+nx*ny]
for i = 2:nsteps
    global truth_full

    truth_full = step(truth_full, three_level_rhs)
    truth[:,i] = truth_full[1:nx+nx*ny]
end

#-------------------------------------------------------------------------------
# %% Extract and perturb observations
#-------------------------------------------------------------------------------

println("Extracting observations...")

# Make observations
y = H(truth)

# Define observational error covariance matrix
R = σ² * Matrix{Float64}(I,obsdim,obsdim)

# Perturb observations
for i in 1:nsteps
    y[:,i] += √(σ²)*randn(size(y[:,i]))
end

#-------------------------------------------------------------------------------
# %% Define ensemble
#-------------------------------------------------------------------------------

println("Generating ensemble...")

ensemble = gen_ensemble(truth)

#-------------------------------------------------------------------------------
# %% Set up output
#-------------------------------------------------------------------------------

nc = set_up_output()

#-------------------------------------------------------------------------------
# %% Run filter
#-------------------------------------------------------------------------------

println("Running filter...")

for i in 1:nsteps
    println("Step $i")

    if i % assim_freq == 0
        ensrf_assimilate!(ensemble, y[:,i], R)
    end

    # Write output
    if i % write_freq == 0
        println("Writing output")
        output(nc, i/write_freq, ensemble, truth[:,i])
    end

    # Forecast step
    for mem in 1:nens
        # Generate stochastic term

        # Step ensemble member forward
        ensemble[:,mem] = step(ensemble[:,mem], two_level_rhs)
    end
end
