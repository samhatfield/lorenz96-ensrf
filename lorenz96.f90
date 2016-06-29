module lorenz96
    use params
    use utils
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
        {% for type in types %}
        module procedure step_param_z_{{ type.name }}
        {% endfor %}
    end interface

    public :: ode_param_z
    interface ode_param_z
        {% for type in types %}
        module procedure ode_param_z_{{ type.name }}
        {% endfor %}
    end interface

    public :: dXdT
    interface dXdT
        {% for type in types %}
        module procedure dXdT_{{ type.name }}
        {% endfor %}
    end interface

    public :: dYdT_param_z
    interface dYdT_param_z
        {% for type in types %}
        module procedure dYdT_param_z_{{ type.name }}
        {% endfor %}
    end interface

    contains
        !===========================================================================
        ! Full three-level model (for full two-level model, just set n_z = 0)
        !===========================================================================

        ! Step forward once
        pure function step(prev_state)
            real(dp), dimension(truth_dim), intent(in) :: prev_state
            real(dp), dimension(truth_dim) :: step, k1, k2, k3, k4

            ! 4th order Runge-Kutta
            k1 = ode(prev_state)
            k2 = ode(prev_state+0.5_dp*dt*k1)
            k3 = ode(prev_state+0.5_dp*dt*k2)
            k4 = ode(prev_state+dt*k3) 

            step = prev_state + (dt/6._dp)*(k1 + 2._dp*k2 + 2._dp*k3 + k4)
        end function
        
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
        end function

        ! X ODE
        {% for type in types %}
        pure function dXdT_{{ type.name }}(x, y) result(dXdT)
            {{ type.code }}, dimension(n_x), intent(in) :: x
            {{ type.code }}, dimension(n_x*n_y), intent(in) :: y
            {{ type.code }}, dimension(n_x) :: dXdT
            {{ type.code }}, dimension(n_x) :: sum_y
            {{ type.code }} :: h_, c_, b_, g_X_, f_

            h_ = h
            c_ = c
            b_ = b
            g_X_ = g_X
            f_ = f

            ! Sum all y's for each x, making an n_x length vector of y sums
            sum_y = sum_2d(reshape(y, (/n_y,n_x/)))

            dXdT = cshift(x,-1)*(cshift(x,1)-cshift(x,-2)) - g_X_*x
            dXdT = dXdT + f_ - (h_*c_/b_)*sum_y
        end function
        {% endfor %}

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
        end function

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
        end function
        
        !===========================================================================
        ! Three-level model with parametrised Z dynamics
        !===========================================================================     

        ! Step forward once
        {% for type in types %}
        function step_param_z_{{ type.name }}(prev_state, stoch) result(step)
            {{ type.code }}, dimension(n_x+n_x*n_y), intent(in) :: prev_state
            {{ type.code }}, dimension(n_x*n_y), intent(in) :: stoch
            {{ type.code }}, dimension(n_x+n_x*n_y) :: step, k1, k2, k3, k4  
            {{ type.code }} :: half, six, two, dt_rpe

            half = 0.5_dp
            six = 6.0_dp
            two = 2.0_dp
            dt_rpe = dt
            
            ! 4th order Runge-Kutta
            k1 = ode_param_z(prev_state, stoch)
            k2 = ode_param_z(prev_state+half*dt_rpe*k1, stoch)
            k3 = ode_param_z(prev_state+half*dt_rpe*k2, stoch)
            k4 = ode_param_z(prev_state+dt_rpe*k3, stoch)

            step = prev_state + (dt_rpe/six)*(k1 + two*k2 + two*k3 + k4)
        end function

        ! The three-level system of ODEs for the Lorenz '96 system, with
        ! parametrized Z
        function ode_param_z_{{ type.name }}(state, stoch) result(ode)
            {{ type.code }}, dimension(state_dim), intent(in) :: state
            {{ type.code }}, dimension(n_x*n_y), intent(in) :: stoch
            {{ type.code }}, dimension(n_x) :: x
            {{ type.code }}, dimension(n_x*n_y) :: y
            {{ type.code }}, dimension(state_dim) :: ode

            ! Break up state vector into components
            x = state(:n_x)
            y = state(n_x+1:)

            ! Find derivative of each component separately
            ode(:n_x) = dXdT(x, y)
            ode(n_x+1:) = dYdT_param_z(x, y, stoch)
        end function
        
        ! Y ODE, parametrized Z
        function dYdT_param_z_{{ type.name }}(x, y, stoch) result(dYdT)
            {{ type.code }}, dimension(n_x), intent(in) :: x
            {{ type.code }}, dimension(n_x*n_y), intent(in) :: y
            {{ type.code }}, dimension(n_x*n_y), intent(in) :: stoch
            {{ type.code }}, dimension(n_x*n_y) :: dYdT
            {{ type.code }}, dimension(n_x*n_y) :: x_rpt
            {{ type.code }}, dimension(n_x*n_y) :: tend_z
            {{ type.code }}, dimension(5) :: coeffs
            {{ type.code }} :: h_, c_, b_, g_Y_
            integer :: k

            coeffs(1) = 0.001892_dp
            coeffs(2) = -0.066811_dp
            coeffs(3) = 0.131826_dp
            coeffs(4) = 0.242742_dp
            coeffs(5) = 0.039970_dp

            h_ = h
            c_ = c
            b_ = b
            g_Y_ = g_Y

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            ! Compute the Z tendency from the chosen parametrization scheme
            ! Deterministic, 4th order polynomial
            tend_z = (coeffs(1)*y**4) + (coeffs(2)*y**3) + &
                & (coeffs(3)*y**2) + (coeffs(4)*y) + coeffs(5)
                
            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - g_Y_*y
            dYdT = dYdT + (h_*c_/b_)*x_rpt        
            dYdT = dYdT - (tend_z + stoch)
        end function
        {% endfor %}
end module lorenz96
