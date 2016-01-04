import numpy as np
from numpy.linalg import inv

RHO = 1.4

def assimilate(ensemble, obs_vec, Ro):
    ens_mean = np.mean(ensemble, axis=0)
    num_ens = len(ensemble)

    ens_trans = np.transpose(ensemble)

    obs_table = np.transpose([obs_vec,] * num_ens)

    Ensfp = ens_trans - np.transpose([ens_mean,] * num_ens)

    HEnsfp = Ensfp

    covfHt = (RHO/(num_ens-1)) * np.dot(Ensfp, np.transpose(HEnsfp))

    HcovfHt = (RHO/(num_ens-1)) * np.dot(HEnsfp, np.transpose(HEnsfp)) + Ro

    Gain = np.dot(covfHt, inv(HcovfHt))

    ens_trans = ens_trans + np.dot(Gain, obs_table - ens_trans)

    return [x for x in np.transpose(ens_trans)]
