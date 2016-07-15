!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> An observing system simulation experiment (OSSE) using the three-tier Lorenz
!> '96/'95 system as the truth, a two tier system as the model, and the
!> ensemble square root filter to assimilate observations. A reduced precision
!> emulator is included, allowing the precision of the forecast and update
!> cycles to be reduced.
!> The token 'PRECISION' is substituted with the specified precision during
!> compilation. The truth and observations are always in double precision.
!>
!> Bibliography
!> ============
!> Thornes et al., "On the Use of Scale-Dependent Precision in Earth-
!> System Modelling", Q. J. Roy. Meteor. Soc. 2016 (in print)
!> Whitaker and Hamill, "Ensemble Data Assimilation without Perturbed
!> Observations", Mon. Weather. Rev. 130 (2002)
program main
    use params
    use lorenz96, only: step, step_param_z
    use utils, only: randn, ar_1, std, identity, real
    use analysis, only: ensrf_assimilate
    use setup, only: write_params, spin_up, gen_ensemble, time_seed
    use rp_emulator
    use observation, only: observe

    implicit none

    !===========================================================================
    ! Declare variables
    !===========================================================================

    integer, parameter :: file_1 = 20

    ! Loop counters
    integer :: i, j

    ! Truth run, observations, ensemble and observation covariance matrix
    real(dp) :: truth(truth_dim, n_steps)
    real(dp) :: obs(obs_dim, n_steps)
    PRECISION :: ensemble(state_dim, n_ens)
    real(dp) :: obs_covar(obs_dim, obs_dim)

    ! For storing norms of each ensemble member (used for output)
    real(dp) :: x_norms(n_ens)

    ! Stores stochastic components for each ensemble member
    PRECISION :: stochs(n_x*n_y, n_ens)

    ! Set up reduced precision emulator
    RPE_DEFAULT_SBITS = sbits
    RPE_IEEE_HALF = .true.

    ! Initialize stochastic term matrix
    stochs(:, :) = 0.0_dp

    ! Seed RNG
    call time_seed()

    !===========================================================================
    ! Spin up
    !===========================================================================

    print *, "Spinning up..."

    truth(:, 1) = spin_up()

    !===========================================================================
    ! Truth run
    !===========================================================================

    print *, "Generating truth..."

    do i = 2, n_steps
        truth(:, i) = step(truth(:, i-1))
    end do

    !===========================================================================
    ! Extract and perturb observations
    !===========================================================================

    print *, "Extracting observations..."

    ! Make observations
    obs = observe(truth) 

    ! Define observational error covariance matrix
    obs_covar = y_var * identity(obs_dim)

    ! Perturb observations
    do i = 1, n_steps
        do j = 1, obs_dim
            obs(j, i) = obs(j, i) + randn(0.0_dp, sqrt(y_var))
        end do
    end do

    !===========================================================================
    ! Define ensemble
    !===========================================================================

    print *, "Generating ensemble..."

    ! Generate an ensemble by taking random samples from truth run
    ! (Equivalent to sampling the climatology)
    ensemble = gen_ensemble(truth)

    !===========================================================================
    ! Write metadata to top of output file
    !===========================================================================

    call write_params()

    !===========================================================================
    ! Run filter
    !===========================================================================

    print *, "Running filter..."

    open(unit=file_1, file="results.yml", action="write", position="append")
    do i = 1, n_steps
       ! Print every 1000th timestep
        if (mod(i, 1000) == 0) then
            write(*,*) 'Step ', i 
        end if

        ! Analysis step (only run every assim_freq steps)
        if (mod(i, assim_freq) == 0) then
            ensemble = ensrf_assimilate(ensemble, obs(:, i), obs_covar)
        end if

        ! Write norm and std of X norms, rms forecast error and truth and
        ! observation vector norm for this timestep
        x_norms = norm2(real(ensemble(:n_x,:)), 1)
        write (file_1, '(6f11.6)') sum(x_norms)/real(n_ens, dp), std(x_norms), &
            & norm2(truth(:n_x,i))

        ! Forecast step
        do j = 1, n_ens
            ! Generate stochastic term
            stochs(:, j) = ar_1(stochs(:, j))

            ! Step ensemble member forward
            ensemble(:, j) = step_param_z(ensemble(:, j), stochs(:, j))
        end do
    end do
    close(file_1)
end program main
