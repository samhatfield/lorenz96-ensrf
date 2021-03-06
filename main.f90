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
    use setup, only: spin_up, gen_ensemble, time_seed
    use rp_emulator
    use observation, only: observe
    use io, only: setup_output, output, open_file, close_file

    implicit none

    !===========================================================================
    ! Declare variables
    !===========================================================================

    ! Loop counters
    integer :: i, j

    ! Truth run, observations, ensemble and observation covariance matrix
	! 'truth_full' stores the full truth state, X, Y and Z variables
	! 'truth' just stores the X and Y variables (at every timestep) to save on memory
    real(dp) :: truth(state_dim, n_steps)
    real(dp) :: truth_full(truth_dim)
    real(dp) :: obs(obs_dim, n_steps)
    PRECISION :: ensemble(state_dim, n_ens)
    real(dp) :: obs_covar(obs_dim, obs_dim)

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

    truth_full = spin_up()
    truth(:,1) = truth_full(:n_x+n_x*n_y)

    !===========================================================================
    ! Truth run
    !===========================================================================

    print *, "Generating truth..."

    do i = 2, n_steps
        truth_full = step(truth_full)
        truth(:,i) = truth_full(:n_x+n_x*n_y)
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
    ! Setup output
    !===========================================================================

    call setup_output()

    !===========================================================================
    ! Run filter
    !===========================================================================

    print *, "Running filter..."

    call open_file()
    do i = 1, n_steps
       ! Print every 1000th timestep
        if (mod(i, 1000) == 0) then
            print *, 'Step ', i 
        end if

        ! Analysis step (only run every assim_freq steps)
        if (mod(i, assim_freq) == 0) then
            ensemble = ensrf_assimilate(ensemble, obs(:, i), obs_covar)
        end if

        ! Write output
        if (mod(i, write_freq) == 0) then
            call output(ensemble, truth(:, i), i / write_freq)
        end if

        ! Forecast step
        do j = 1, n_ens
            ! Generate stochastic term
            stochs(:, j) = ar_1(stochs(:, j))

            ! Step ensemble member forward
            ensemble(:, j) = step_param_z(ensemble(:, j), stochs(:, j))
        end do
    end do
    call close_file()
end program main
