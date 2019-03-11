include("params.jl")
include("lorenz96.jl")

using Random: randn

"Spins up the truth run."
function spinup()
    # Some random initial conditions (it doesn't really matter)
    truth = 8.0*ones(truthdim)
    truth[nx+1:nx+nx*ny] = 0.5*randn(nx*ny)
    truth[nx+nx*ny+1:end] = 0.5*randn(nx*ny*nz)
    truth[4] = 8.008

    # Spin up
    for i = 1:5000
        truth = step(truth, three_level_rhs)
    end

    return truth
end

"Generate the first background ensemble by sampling from the truth run."
function gen_ensemble(truth)
    ensemble = zeros(statedim, nens)

    for mem in 1:nens
        ensemble[:,mem] = truth[:,ceil(Int, nsteps*rand())]
    end

    return ensemble
end
