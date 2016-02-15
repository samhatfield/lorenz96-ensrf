program main
    use params
    use lorenz96, only: step
    use utils, only: randn, time_seed
    use analysis
    use metadata

    implicit none

    !===========================================================================
    ! Declare globals
    !===========================================================================

    integer, parameter :: file_1 = 20
    integer, parameter :: file_2 = 21
    integer, parameter :: file_3 = 22

    ! Loop counters
    integer :: i, j

    real(dp), dimension(state_dim) :: initial_truth
    real(dp), dimension(state_dim, n_steps) :: truth_run
    real(dp), dimension(obs_dim, n_steps) :: observations
    real(dp), dimension(state_dim, n_ens) :: ensemble
    real(dp), dimension(obs_dim, obs_dim) :: obs_covar

    ! For storing norms of each ensemble member (used for output)
    real(dp), dimension(n_ens) :: x_norms

    ! For storing RMS forecast error, and ensemble mean vector
    real(dp) :: rms_err
    real(dp), dimension(state_dim) :: ens_mean

    ! Seed
    call time_seed()

    !===========================================================================
    ! Spin up
    !===========================================================================

    write(*,*) "Spinning up..."

    ! Initial conditions for spin up
    initial_truth(:n_x) = (/ (8, i = 1, n_x) /)
    initial_truth(n_x+1:n_x+n_x*n_y) = (/ (randn(0d0, 0.5d0), i = 1, n_x*n_y) /)
    initial_truth(n_x+n_x*n_y+1:) = (/ (randn(0d0, 0.5d0), i = 1, n_x*n_y*n_z) /)
    initial_truth(4) = 8.008

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
    observations = observe(truth_run) 

    ! Define observational error covariance matrix (diagonal matrix of variances)
    forall(i = 1:obs_dim, j = 1:obs_dim) obs_covar(i, j) = var_obs * (i/j)*(j/i)

    ! Perturb observations
    do i = 1, n_steps
        do j = 1, obs_dim
            ! Observation perturbation has variance of ~0.1
            observations(j, i) = observations(j, i) + randn(0.0d0, sig_obs)
        end do
    end do

    !===========================================================================
    ! Define ensemble
    !===========================================================================

    write(*,*) "Generating ensemble..."

    ! Perturb initial truth to generate members
    do i = 1, n_ens
        do j = 1, state_dim
            ! Member perturbation has variance of ~3
            ensemble(j, i) = initial_truth(j) + randn(0.0d0, 1.73d0)
        end do
    end do

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
        ! Print every 100th timestep
        if (mod(i, 100) == 0) then
            write(*,*) 'Step ', i 
        end if

        ! Forecast step
        do j = 1, n_ens
            ensemble(:, j) = step(ensemble(:, j))
        end do

        ! Analysis step
        if (mod(i, assim_freq) == 0) then
            ensemble = assimilate(ensemble, observations(:, i), obs_covar, sig_obs)
        end if

        ! Write upper, lower and average norm of ensemble members, rms forecast
        ! error and truth and observation vector norm for this timestep
        x_norms = norm2(ensemble(:n_x,:), 1)
        ens_mean = (/ (sum(ensemble(i, :))/real(n_ens) , i = 1, state_dim) /)
        rms_err = norm2(truth_run(:n_x, i) - ens_mean(:n_x))

        write (file_2, '(6f11.6)') maxval(x_norms), sum(x_norms)/real(n_ens), &
            & minval(x_norms), norm2(truth_run(:n_x,i)), norm2(observations(:,i)), &
            & rms_err
    end do
    close(file_2)
end program main
