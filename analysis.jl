include("params.jl")
include("observation.jl")

using Statistics: mean

#= Performs sequential ensemble square-root filtering data assimilation with the given background
   ensemble, observation vector and observation covariance matrix. =#
function ensrf_assimilate!(Pf, y, R)
    # Mean ensemble vector
    Pf̅ = dropdims(mean(Pf, dims=2), dims=2)

    # Form the background ensemble perturbation matrix (with covariance inflation)
    Xf = ρ*(Pf .- Pf̅)

    for i in 1:obsdim
        # Ensemble covariance times transpose of observation matrix
        PfHᵀ = Xf * H(Xf, i) / (nens-1)

        HPfHᵀ = H(PfHᵀ,i)

        # Kalman gain
        K = PfHᵀ / (HPfHᵀ + R[i,i])

        # Localization
        # if (loc >= 0.0)
        #
        # end

        # Update ensemble mean
        Pf̅ += K*(y[i] - H(Pf̅,i))

        # Update perturbation
        α = 1.0/(1.0 + √(R[i,i]/HPfHᵀ + R[i,i]))
        for j in 1:nens
            Xf[:,j] -= α*K*H(Xf[:,j],i)
        end
    end

    # Form final ensemble
    for i in 1:nens
        Pf[:,i] = Pf̅ + Xf[:,i]
    end
end