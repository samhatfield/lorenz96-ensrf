program main
    use params
    use lorenz96, only: step, step_param_z
    use utils, only: randn, time_seed, additive_noise
    use analysis
    use metadata
    use rp_emulator

    implicit none

    !===========================================================================
    ! Declare globals
    !===========================================================================

    integer, parameter :: file_1 = 20
    integer, parameter :: file_2 = 21
    integer, parameter :: file_3 = 22

    ! Loop counters
    integer :: i, j

    real(dp), dimension(truth_dim) :: initial_truth
    real(dp), dimension(truth_dim, n_steps) :: truth_run
    real(dp), dimension(obs_dim, n_steps) :: obs
    real(dp), dimension(state_dim, n_ens) :: ensemble
    real(dp), dimension(obs_dim, obs_dim) :: obs_covar
    real(dp), dimension(truth_dim) :: climatology_mean
    real(dp), dimension(truth_dim) :: climatology_std

    ! For storing norms of each ensemble member (used for output)
    real(dp), dimension(n_ens) :: x_norms

    ! For storing RMS forecast error, and ensemble mean vector
    real(dp) :: rms_err
    real(dp), dimension(state_dim) :: ens_mean
    
    ! Stores stochastic components for each ensemble members
    real(dp), dimension(n_x*n_y, n_ens) :: stochs

    ! Seed
    call time_seed()

    !===========================================================================
    ! Spin up
    !===========================================================================

    write(*,*) "Spinning up..."

    ! Initial conditions for spin up
    initial_truth(:n_x) = (/ (8, i = 1, n_x) /)
    initial_truth(n_x+1:n_x+n_x*n_y) = (/ (randn(0._dp, 0.5_dp), i = 1, n_x*n_y) /)
    initial_truth(n_x+n_x*n_y+1:) = (/ (randn(0._dp, 0.5_dp), i = 1, n_x*n_y*n_z) /)
    initial_truth(4) = 8.008_dp

    ! Spin up
    do i = 1, 5000
        initial_truth = step(initial_truth)
    end do

    !===========================================================================
    ! Truth run
    !===========================================================================

    write(*,*) "Generating truth..."

    truth_run(:, 1) = initial_truth
    do i = 2, n_steps
        truth_run(:, i) = step(truth_run(:, i-1))
    end do

    !===========================================================================
    ! Extract and perturb observations
    !===========================================================================

    write(*,*) "Extracting observations..."

    ! Make observations
    obs = observe(truth_run) 

    ! Define observational error covariance matrix (diagonal matrix of variances)
    forall(i = 1:obs_dim, j = 1:obs_dim) obs_covar(i, j) = y_var * (i/j)*(j/i)

    ! Perturb observations
    do i = 1, n_steps
        do j = 1, obs_dim
            obs(j, i) = obs(j, i) + randn(0.0_dp, sqrt(y_var))
        end do
    end do

    !===========================================================================
    ! Define ensemble
    !===========================================================================

    write(*,*) "Generating ensemble..."

    ! Get climatology from time average of truth
    climatology_mean = (/ (sum(truth_run(j, :))/real(n_steps, dp), j = 1, truth_dim) /)
    climatology_std = (/ (std(truth_run(j, :)), j = 1, truth_dim) /)

    ! Perturb climatology to generate members
    do i = 1, n_ens
        do j = 1, state_dim
            ensemble(j, i) = climatology_mean(j) + randn(0.0_dp, climatology_std(j))
        end do
    end do

    !===========================================================================
    ! Initialise stochastic components vector
    !===========================================================================

    stochs(:, :) = 0.0_dp

    !===========================================================================
    ! Write metadata to top of output file
    !===========================================================================

    call write_params()

    !===========================================================================
    ! Run filter
    !===========================================================================

    write(*,*) "Running filter..."

    open(unit=file_2, file="results.yml", action="write", position="append")
    do i = 1, n_steps
        ! Write upper, lower and average norm of ensemble members, rms forecast
        ! error and truth and observation vector norm for this timestep
        x_norms = norm2(ensemble(:n_x,:), 1)
        ens_mean = (/ (sum(ensemble(j, :))/real(n_ens, dp) , j = 1, state_dim) /)
        rms_err = norm2(truth_run(:n_x, i) - ens_mean(:n_x))

        write (file_2, '(6f11.6)') sum(x_norms)/real(n_ens, dp), std(x_norms), &
            & norm2(truth_run(:n_x,i)), norm2(obs(:,i)), rms_err
    
        ! Print every 100th timestep
        if (mod(i, 100) == 0) then
            write(*,*) 'Step ', i 
        end if

        ! Forecast step
        do j = 1, n_ens
            stochs(:, j) = additive_noise(stochs(:, j))
            ensemble(:, j) = step_param_z(ensemble(:, j), stochs(:, j))
        end do

        ! Analysis step
        if (mod(i, assim_freq) == 0) then
            ensemble = assimilate(ensemble, obs(:, i), obs_covar)
        end if
    end do
    close(file_2)
end program main
