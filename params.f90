module params
    implicit none
    
    integer, parameter :: dp = kind(0.d0)
    integer, parameter :: sp = kind(0.0)

    !Start and end times, and number of time steps
    real(dp), parameter :: dt = 0.005_dp
    real(dp), parameter :: fin = 12_dp
    integer, parameter :: n_steps = fin / dt

    ! State vector properties
    integer, parameter :: n_x = 8
    integer, parameter :: n_y = 8
    integer, parameter :: n_z = 8
    
    ! State vector dimension of full model
    integer, parameter :: truth_dim = n_x + n_x*n_y + n_x*n_y*n_z
    
    ! State vector dimension of parametrised model
    integer, parameter :: state_dim = n_x + n_x*n_y
    
    ! Number of ensemble members
    integer, parameter :: n_ens = 40

    ! Frequency of assimilations, i.e. 1 = every timestep
    integer, parameter :: assim_freq = 1

    ! Observation error variance
    real(dp), parameter :: y_var = 0.1_dp

    ! Dimension of observation vector
    integer, parameter :: obs_dim = n_x*n_y

    ! Default significand
    integer, parameter :: sbits = 10

    ! Covariance inflation factor
    real(dp), parameter :: rho = 1.003_dp

    ! Covariance localisation length scales
    real(dp), parameter :: loc = 1.2_dp
    real(dp), parameter :: loc_y = loc*26
    real(dp), parameter :: loc_x = loc*8
end module params
