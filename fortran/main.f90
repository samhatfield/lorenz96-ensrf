program main
    implicit none

    !===========================================================================
    ! Declare globals
    !===========================================================================

    ! Start and end times, and number of time steps
    real, parameter :: dt = 0.005
    real, parameter :: fin = 10
    integer, parameter :: n_steps = fin / dt

    ! Loop counters
    integer :: i, j

    ! Number of ensemble members
    integer, parameter :: n_ens = 500

    ! Frequency of assimilations, i.e. 1 = every timestep
    integer, parameter :: assim_freq = 1

    ! Model parameters
    integer, parameter :: n_x = 8
    integer, parameter :: n_y = 32
    integer, parameter :: state_dim = n_x + n_x * n_y
    real, parameter :: f = 20
    real :: h = 1
    real :: c = 4
    real :: b = 10

    ! Observation error variance
    real, parameter :: var_obs = 0.1
    real, parameter :: sigma_obs = sqrt(var_obs)

    ! For now, hard code the observation state vector dimension, assuming that
    ! we are only assimilating the Y variables
    ! Later on, I'll find a way for it to automatically calculate the
    ! observation state vector dimension at compile time
    integer, parameter :: obs_dim = n_y

    real, dimension(state_dim) :: initial_truth
    real, dimension(state_dim, n_steps) :: truth_run
    real, dimension(obs_dim, n_steps) :: observations
    real, dimension(state_dim, n_ens) :: ensemble
    real, dimension(obs_dim, obs_dim) :: obs_covar

    !===========================================================================
    ! Spin up
    !===========================================================================

    !write(*,*) "Spinning up..."

    ! Initial conditions for spin up
    initial_truth(:n_x) = (/ (8, i = 1, n_x) /)
    initial_truth(n_x+1:) = (/ (0.5, i = 1, n_x*n_y) /)
    initial_truth(4) = 8.008

    ! Spin up
    do i = 1, 5000
        initial_truth = step(initial_truth)
    end do

    !===========================================================================
    ! Truth run
    !===========================================================================

    !write(*,*) "Generating truth..."

    truth_run(:, 1) = initial_truth
    do i = 2, n_steps
        truth_run(:, i) = step(truth_run(:, i-1))
    end do

    !===========================================================================
    ! Extract and perturb observations
    !===========================================================================

    !write(*,*) "Extracting observations..."

    ! Make observations
    observations = observe(truth_run) 

    ! Define observational error covariance matrix (diagonal matrix of variances)
    forall(i = 1:obs_dim, j = 1:obs_dim) obs_covar(i, j) = var_obs * (i/j)*(j/i)

    ! Perturb observations
    do i = 1, n_steps
        do j = 1, obs_dim
            ! Member perturbation has variance of ~3
            observations(j, i) = observations(j, i) + randn(0.0, sigma_obs)
        end do
    end do

    !===========================================================================
    ! Define ensemble
    !===========================================================================

    !write(*,*) "Generating ensemble..."

    ! Perturb initial truth to generate members
    do i = 1, n_ens
        do j = 1, state_dim
            ! Member perturbation has variance of ~3
            ensemble(j, i) = initial_truth(j) + randn(0.0, 1.73)
        end do
    end do

    !===========================================================================
    ! Run filter
    !===========================================================================

    !write(*,*) "Running filter..."

    do i = 1, n_steps
        ! Print every 100th timestep
        if (mod(i, 100) == 0) then
            !write(*,*) 'Step ', i 
        end if

        ! Forecast step
        do j = 1, n_ens
            ensemble(:, j) = step(ensemble(:, j))
        end do

!        ! Analysis step
!        if (mod(i, assim_freq) == 0) then
!            ensemble = assimilate(ensemble, observations(i), Ro, sigma_o)
!        end if

        print "(264F10.5)", (/ (sum(ensemble(i, :))/real(n_ens) , i = 1, state_dim) /)
    end do

    contains
        ! Step forward once
        pure function step(prev_state)
            real, dimension(state_dim), intent(in) :: prev_state
            real, dimension(n_x) :: x, k1, k2, k3, k4
            real, dimension(n_x*n_y) :: y, l1, l2, l3, l4
            real, dimension(state_dim) :: step

            x = prev_state(:n_x)
            y = prev_state(n_x+1:)
        
            ! 2 coupled ODEs 4th order Runge-Kutta
            k1 = dXdT(x, y)
            l1 = dYdT(x, y)
            k2 = dXdT(x+0.5*dt*k1, y+0.5*dt*l1)
            l2 = dYdT(x+0.5*dt*k1, y+0.5*dt*l1)
            k3 = dXdT(x+0.5*dt*k2, y+0.5*dt*l2)
            l3 = dYdT(x+0.5*dt*k2, y+0.5*dt*l2)
            k4 = dXdT(x+dt*k3, y+dt*l3)
            l4 = dYdT(x+dt*k3, y+dt*l3)
        
            x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
            y = y + (dt/6)*(l1 + 2*l2 + 2*l3 + l4)

            ! Concatenate X's and Y's into return value
            step(:n_x) = x
            step(n_x+1:) = y
        end function step

        ! X ODE
        pure function dXdT(x, y)
            real, dimension(n_x), intent(in) :: x
            real, dimension(n_x*n_y), intent(in) :: y
            real, dimension(n_x) :: dXdT
            real, dimension(n_x) :: sum_y

            sum_y = sum(reshape(y, (/n_y,n_x/)), dim=1)

            dXdT = cshift(x,-1)*(cshift(x,1)-cshift(x,-2)) - x + f
            dXdT = dXdT - (h*c/b)*sum_y
        end function dXdT

        ! Y ODE
        pure function dYdT(x, y)
            real, dimension(n_x), intent(in) :: x
            real, dimension(n_x*n_y), intent(in) :: y
            real, dimension(n_x*n_y) :: dYdT
            real, dimension(n_x*n_y) :: x_rpt
            integer :: k

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - c*y
            dYdT = dYdT + (h*c/b)*x_rpt
        end function dYdT

        pure function observe(state)
            real, dimension(:, :), intent(in) :: state
            real, dimension(obs_dim, size(state, 2)) :: observe

            observe = state(n_x+1:,:)
        end function observe

        ! Generates a random number drawn for the specified normal distribution
        function randn(mean, stdev)
            real, intent(in) :: mean, stdev
            real :: randn, rand(2), u, v

            call random_number(rand)

            ! Box-Muller method
            u = (-2.0d0 * log(rand(1))) ** 0.5
            v =   2.0d0 * 6.28318530718 * rand(2)
            randn = mean + stdev * u * sin(v)
        end function randn
end program main
