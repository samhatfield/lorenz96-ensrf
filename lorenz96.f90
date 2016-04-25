module lorenz96
    use params
    use utils, only: randn, sum_2d_rpe
    use rp_emulator

    implicit none

    ! Model parameters
    real(dp), parameter :: f = 20._dp
    real(dp), parameter :: h = 1._dp
    real(dp), parameter :: c = 10._dp
    real(dp), parameter :: b = 10._dp
    real(dp), parameter :: e = 10._dp
    real(dp), parameter :: d = 10._dp
    real(dp), parameter :: g_X = 1._dp
    real(dp), parameter :: g_Y = c
    real(dp), parameter :: g_Z = e

    !===========================================================================
    ! Interfaces for overloaded definitions
    !===========================================================================

    public :: step_param_z
    interface step_param_z
        module procedure step_param_z
        module procedure step_param_z_rpe
    end interface step_param_z

    public :: ode_param_z
    interface ode_param_z
        module procedure ode_param_z
        module procedure ode_param_z_rpe
    end interface ode_param_z

    public :: dXdT
    interface dXdT
        module procedure dXdT
        module procedure dXdT_rpe
    end interface dXdT

    public :: dYdT_param_z
    interface dYdT_param_z
        module procedure dYdT_param_z
        module procedure dYdT_param_z_rpe
    end interface dYdT_param_z

    contains
        !===========================================================================
        ! Full three-level model (for full two-level model, just set n_z = 0)
        !===========================================================================

        ! Step forward once
        function step(prev_state)
            real(dp), dimension(truth_dim), intent(in) :: prev_state
            real(dp), dimension(truth_dim) :: step, k1, k2, k3, k4

            ! 4th order Runge-Kutta
            k1 = ode(prev_state)
            k2 = ode(prev_state+0.5_dp*dt*k1)
            k3 = ode(prev_state+0.5_dp*dt*k2)
            k4 = ode(prev_state+dt*k3) 

            step = prev_state + (dt/6._dp)*(k1 + 2._dp*k2 + 2._dp*k3 + k4)
        end function step
        
        ! The full three-level system of ODEs for the Lorenz '96 system
        pure function ode(state)
            real(dp), dimension(truth_dim), intent(in) :: state
            real(dp), dimension(n_x) :: x
            real(dp), dimension(n_x*n_y) :: y
            real(dp), dimension(n_x*n_y*n_z) :: z
            real(dp), dimension(truth_dim) :: ode

            ! Break up state vector into components
            x = state(:n_x)
            y = state(n_x+1:n_x+n_x*n_y)
            z = state(n_x+n_x*n_y+1:)

            ! Find derivative of each component separately
            ode(:n_x) = dXdT(x, y)
            ode(n_x+1:n_x+n_x*n_y) = dYdT(x, y, z)
            ode(n_x+n_x*n_y+1:) = dZdT(y, z)
        end function ode

        ! X ODE
        pure function dXdT(x, y)
            real(dp), dimension(n_x), intent(in) :: x
            real(dp), dimension(n_x*n_y), intent(in) :: y
            real(dp), dimension(n_x) :: dXdT
            real(dp), dimension(n_x) :: sum_y

            ! Sum all y's for each x, making an n_x length vector of y sums
            sum_y = sum(reshape(y, (/n_y,n_x/)), dim=1)

            dXdT = cshift(x,-1)*(cshift(x,1)-cshift(x,-2)) - g_X*x + f
            dXdT = dXdT - (h*c/b)*sum_y
        end function dXdT

        ! Y ODE
        pure function dYdT(x, y, z)
            real(dp), dimension(n_x), intent(in) :: x
            real(dp), dimension(n_x*n_y), intent(in) :: y
            real(dp), dimension(n_x*n_y*n_z), intent(in) :: z
            real(dp), dimension(n_x*n_y) :: dYdT
            real(dp), dimension(n_x*n_y) :: x_rpt
            real(dp), dimension(n_x*n_y) :: sum_z
            integer :: k

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            ! Sum all z's for each y, making an n_x*n_y length vector of z sums
            sum_z = sum(reshape(z, (/n_z,n_x*n_y/)), dim=1)

            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - g_Y*y
            dYdT = dYdT + (h*c/b)*x_rpt
            dYdT = dYdT - (h*e/d)*sum_z
        end function dYdT

        ! Z ODE
        pure function dZdT(y, z)
            real(dp), dimension(n_x*n_y), intent(in) :: y
            real(dp), dimension(n_x*n_y*n_z), intent(in) :: z
            real(dp), dimension(n_x*n_y*n_z) :: dZdT
            real(dp), dimension(n_x*n_y*n_z) :: y_rpt
            integer :: k

            ! Repeat elements of y n_z times
            y_rpt = (/ (y(1+(k-1)/n_z), k = 1, n_x*n_y*n_z) /)

            dZdT = e*b*cshift(z,-1)*(cshift(z,1)-cshift(z,-2)) - g_Z*Z
            dZdT = dZdT + (h*e/d)*y_rpt
        end function dZdT
        
        !===========================================================================
        ! Three-level model with parametrised Z dynamics
        !===========================================================================     

        ! Step forward once
        pure function step_param_z(prev_state, stoch) result(step)
            real(dp), dimension(n_x+n_x*n_y), intent(in) :: prev_state
            real(dp), dimension(n_x*n_y), intent(in) :: stoch
            real(dp), dimension(n_x+n_x*n_y) :: step, k1, k2, k3, k4  
            
            ! 4th order Runge-Kutta
            k1 = ode_param_z(prev_state, stoch)
            k2 = ode_param_z(prev_state+0.5_dp*dt*k1, stoch)
            k3 = ode_param_z(prev_state+0.5_dp*dt*k2, stoch)
            k4 = ode_param_z(prev_state+dt*k3, stoch)

            step = prev_state + (dt/6._dp)*(k1 + 2._dp*k2 + 2._dp*k3 + k4)
        end function step_param_z

        ! The three-level system of ODEs for the Lorenz '96 system, with
        ! parametrized Z
        pure function ode_param_z(state, stoch) result(ode)
            real(dp), dimension(state_dim), intent(in) :: state
            real(dp), dimension(n_x*n_y), intent(in) :: stoch
            real(dp), dimension(n_x) :: x
            real(dp), dimension(n_x*n_y) :: y
            real(dp), dimension(state_dim) :: ode

            ! Break up state vector into components
            x = state(:n_x)
            y = state(n_x+1:)

            ! Find derivative of each component separately
            ode(:n_x) = dXdT(x, y)
            ode(n_x+1:) = dYdT_param_z(x, y, stoch)
        end function ode_param_z

        ! Y ODE, parametrized Z
        pure function dYdT_param_z(x, y, stoch) result(dYdT)
            real(dp), dimension(n_x), intent(in) :: x
            real(dp), dimension(n_x*n_y), intent(in) :: y
            real(dp), dimension(n_x*n_y), intent(in) :: stoch
            real(dp), dimension(n_x*n_y) :: dYdT
            real(dp), dimension(n_x*n_y) :: x_rpt
            real(dp), dimension(n_x*n_y) :: tend_z
            integer :: k

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            ! Compute the Z tendency from the chosen parametrization scheme
            ! Deterministic, 4th order polynomial
            tend_z = (0.001892_dp*y**4) + (-0.066811_dp*y**3) + &
                & (0.131826_dp*y**2) + (0.242742_dp*y) + 0.039970_dp
                
            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - g_Y*y
            dYdT = dYdT + (h*c/b)*x_rpt        
            dYdT = dYdT - (tend_z + stoch)
        end function dYdT_param_z
        
        !===========================================================================
        ! Three-level model with parametrised Z dynamics and reduced precision
        !===========================================================================     

        ! Step forward once
        function step_param_z_rpe(prev_state, stoch) result(step)
            type(rpe_var), dimension(n_x+n_x*n_y), intent(in) :: prev_state
            type(rpe_var), dimension(n_x*n_y), intent(in) :: stoch
            type(rpe_var), dimension(n_x+n_x*n_y) :: step, k1, k2, k3, k4  
            
            ! 4th order Runge-Kutta
            k1 = ode_param_z(prev_state, stoch)
            k2 = ode_param_z(prev_state+0.5_dp*dt*k1, stoch)
            k3 = ode_param_z(prev_state+0.5_dp*dt*k2, stoch)
            k4 = ode_param_z(prev_state+dt*k3, stoch)

            step = prev_state + (dt/6._dp)*(k1 + 2._dp*k2 + 2._dp*k3 + k4)
        end function step_param_z_rpe

        ! The three-level system of ODEs for the Lorenz '96 system, with
        ! parametrized Z
        pure function ode_param_z_rpe(state, stoch) result(ode)
            type(rpe_var), dimension(state_dim), intent(in) :: state
            type(rpe_var), dimension(n_x*n_y), intent(in) :: stoch
            type(rpe_var), dimension(n_x) :: x
            type(rpe_var), dimension(n_x*n_y) :: y
            type(rpe_var), dimension(state_dim) :: ode

            ! Break up state vector into components
            x = state(:n_x)
            y = state(n_x+1:)

            ! Find derivative of each component separately
            ode(:n_x) = dXdT(x, y)
            ode(n_x+1:) = dYdT_param_z(x, y, stoch)
        end function ode_param_z_rpe
        
        ! X ODE
        pure function dXdT_rpe(x, y) result(dXdT)
            type(rpe_var), dimension(n_x), intent(in) :: x
            type(rpe_var), dimension(n_x*n_y), intent(in) :: y
            type(rpe_var), dimension(n_x) :: dXdT
            type(rpe_var), dimension(n_x) :: sum_y

            ! Sum all y's for each x, making an n_x length vector of y sums
            sum_y = sum_2d_rpe(reshape(y, (/n_y,n_x/)))

            dXdT = cshift(x,-1)*(cshift(x,1)-cshift(x,-2)) - g_X*x + f
            dXdT = dXdT - (h*c/b)*sum_y
        end function dXdT_rpe
        
        ! Y ODE, parametrized Z
        pure function dYdT_param_z_rpe(x, y, stoch) result(dYdT)
            type(rpe_var), dimension(n_x), intent(in) :: x
            type(rpe_var), dimension(n_x*n_y), intent(in) :: y
            type(rpe_var), dimension(n_x*n_y), intent(in) :: stoch
            type(rpe_var), dimension(n_x*n_y) :: dYdT
            type(rpe_var), dimension(n_x*n_y) :: x_rpt
            type(rpe_var), dimension(n_x*n_y) :: tend_z
            integer :: k

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            ! Compute the Z tendency from the chosen parametrization scheme
            ! Deterministic, 4th order polynomial
            tend_z = (0.001892_dp*y**4) + (-0.066811_dp*y**3) + &
                & (0.131826_dp*y**2) + (0.242742_dp*y) + 0.039970_dp
                
            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - g_Y*y
            dYdT = dYdT + (h*c/b)*x_rpt        
            dYdT = dYdT - (tend_z + stoch)
        end function dYdT_param_z_rpe
end module lorenz96
