import numpy as np
from numpy.random import randn
from numpy.linalg import inv

from observation import observe

RHO = 1.4

def assimilate(ensemble, obs_vec, Ro, sigma_o):
    ens_mean = np.mean(ensemble, axis=1)
    num_ens = ensemble.shape[1]
    num_obs = obs_vec.shape[0]

    # Table of perturbed observations. Each column contains the 
    # observation vector for this timestep with a bit of Gaussian-drawn
    # noise (different noise for each member)
    obs_table = np.empty((num_obs, num_ens))
    for i in range(num_ens):
        obs_table[:,i] = obs_vec + sigma_o * randn(num_obs)

    # 'Anomaly table'. Each column is that member's difference from the mean
    Ensfp = ensemble - ens_mean[:,np.newaxis]

    HEnsfp = observe(Ensfp)

    covfHt = (RHO/(num_ens-1)) * np.dot(Ensfp, np.transpose(HEnsfp))

    HcovfHt = (RHO/(num_ens-1)) * np.dot(HEnsfp, np.transpose(HEnsfp)) + Ro

    Gain = np.dot(covfHt, inv(HcovfHt))

    return ensemble + np.dot(Gain, obs_table - observe(ensemble))
