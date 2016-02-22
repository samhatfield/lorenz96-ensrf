module params
    implicit none
    
    integer, parameter :: dp = kind(1.d0)

    !Start and end times, and number of time steps
    real(dp), parameter :: dt = 0.002
    real(dp), parameter :: fin = 12
    integer, parameter :: n_steps = fin / dt

    ! State vector properties
    integer, parameter :: n_x = 8
    integer, parameter :: n_y = 8
    integer, parameter :: n_z = 8
    integer, parameter :: state_dim = n_x + n_x*n_y + n_x*n_y*n_z

    ! Number of ensemble members
    integer, parameter :: n_ens = 500

    ! Frequency of assimilations, i.e. 1 = every timestep
    integer, parameter :: assim_freq = 1

    ! Observation error variance
    real(dp), parameter :: var_obs = 0.1
    real(dp), parameter :: sig_obs = sqrt(var_obs)

    ! Inflation factor for tuning assimilation (see analysis.f90)
    real(dp), parameter :: rho = 1.4

    ! For now, hard code the observation state vector dimension, assuming that
    ! we are only assimilating the Z variables
    ! Later on, I'll find a way for it to automatically calculate the
    ! observation state vector dimension at compile time
    ! TEMPORARILY OBSERVING ENTIRE STATE
    integer, parameter :: obs_dim = state_dim
end module params
