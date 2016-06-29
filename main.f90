program main
    use params
    use lorenz96, only: step, step_param_z
    use utils, only: randn, time_seed, ar_1
    use analysis
    use metadata
    use rp_emulator
    use observation, only: observe

    implicit none

    !===========================================================================
    ! Declare globals
    !===========================================================================

    integer, parameter :: file_1 = 20

    ! Loop counters
    integer :: i, j

    real(dp), dimension(truth_dim) :: initial_truth
    real(dp), dimension(truth_dim, n_steps) :: truth_run
    real(dp), dimension(obs_dim, n_steps) :: obs
    PRECISION, dimension(state_dim, n_ens) :: ensemble
    real(dp), dimension(obs_dim, obs_dim) :: obs_covar
    PRECISION, dimension(truth_dim) :: climatology_mean
    PRECISION, dimension(truth_dim) :: climatology_std
    real(dp) :: rand

    ! For storing norms of each ensemble member (used for output)
    real(dp), dimension(n_ens) :: x_norms

    ! Stores stochastic components for each ensemble member
    PRECISION, dimension(n_x*n_y, n_ens) :: stochs

    ! Literals
    PRECISION :: zero
    zero = 0.0_dp

    RPE_DEFAULT_SBITS = sbits

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

    ! Generate an ensemble by taking random samples from truth run
    ! (Equivalent to sampling the climatology)
    do i = 1, n_ens
        call random_number(rand)
        j = ceiling(rand * n_steps)
        ensemble(:, i) = truth_run(:n_x+n_x*n_y, j)
    end do

    !===========================================================================
    ! Initialise stochastic components vector
    !===========================================================================

    stochs(:, :) = zero

    !===========================================================================
    ! Write metadata to top of output file
    !===========================================================================

    call write_params()

    !===========================================================================
    ! Run filter
    !===========================================================================

    write(*,*) "Running filter..."

    open(unit=file_1, file="results.yml", action="write", position="append")
    do i = 1, n_steps
       ! Print every 100th timestep
        if (mod(i, 100) == 0) then
            write(*,*) 'Step ', i 
        end if

        ! Analysis step
        if (mod(i, assim_freq) == 0) then
            ensemble = ensrf_assimilate(ensemble, obs(:, i), obs_covar)
        end if

        ! Write norm and std of X norms, rms forecast error and truth and
        ! observation vector norm for this timestep
        x_norms = norm2(real(ensemble(:n_x,:)), 1)

        write (file_1, '(6f11.6)') sum(x_norms)/real(n_ens, dp), std(x_norms), &
            & norm2(truth_run(:n_x,i))

        ! Forecast step
        do j = 1, n_ens
            stochs(:, j) = ar_1(stochs(:, j))
            ensemble(:, j) = step_param_z(ensemble(:, j), stochs(:, j))
        end do
    end do
    close(file_1)
end program main
