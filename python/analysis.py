import numpy as np
from numpy import dot
from numpy.random import randn
from numpy.linalg import inv
from math import sqrt

from observation import observe as H

RHO = 1

def enkf_assimilate(ensemble, obs_vec, Ro, sigma_o):
    # Mean ensemble vector
    ens_mean = np.mean(ensemble, axis=1)

    # Number of ensemble members
    N = ensemble.shape[1]

    # Dimension of observed vector
    num_obs = obs_vec.shape[0]

    # Table of perturbed observations. Each column contains the 
    # observation vector for this timestep with a bit of Gaussian-drawn
    # noise (different noise for each member)
    obs_table = np.empty((num_obs, N))
    for i in range(N):
        obs_table[:,i] = obs_vec + sigma_o * randn(num_obs)

    # 'Anomaly table'. Each column is that member's difference from the mean
    anom_table = ensemble - ens_mean[:,np.newaxis]

    # Ensemble covariance times transpose of observation matrix
    ens_cov_H_T = (RHO/(N-1)) * dot(anom_table, np.transpose(H(anom_table)))

    # Kalman gain
    gain = dot(ens_cov_H_T, inv(H(ens_cov_H_T) + Ro))

    # Combine previous ensemble and observations with Kalman gain to generate
    # next ensemble
    return ensemble + dot(gain, obs_table - H(ensemble))

def ensrf_assimilate(ensemble, obs_vec, Ro, sigma_o):
    # Mean ensemble vector
    ens_mean = np.mean(ensemble, axis=1)

    # Number of ensemble members
    m = ensemble.shape[1]

    # Observation dimension
    p = obs_vec.shape[0]

    # Forecast perturbation matrix
    Z_f = ensemble - ens_mean[:,np.newaxis]

    # Ensemble covariance times transpose of observation matrix
    ens_cov_H_T = dot(Z_f, np.transpose(H(Z_f))) / (m-1)

    # Kalman gain
    gain = dot(ens_cov_H_T, inv(H(ens_cov_H_T) + Ro))

    # Compute analysis mean
    ens_mean += dot(gain, obs_vec - H(ens_mean))

    V = H(Z_f)

    # Assimilate observations serially, as they are uncorrelated
    for i in range(p):
        v = V[i,:]
        v_dot = dot(v, v)

        D = v_dot + sigma_o[i]**2

        beta = 1.0/(D + sqrt(D*sigma_o[i]**2))

        Z_f = (1 - beta*v_dot) * Z_f


    return ens_mean[:,np.newaxis] + Z_f
