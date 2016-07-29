!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions for integrating the Lorenz '96/'95 system with three
!> levels and two-levels (with a parametrized third layer). Jinja2 templates
!> are used to generate overloaded functions of different precisions.
!>
!> Bibliography
!> ============
!> Thornes et al., "On the Use of Scale-Dependent Precision in Earth-
!> System Modelling", Q. J. Roy. Meteor. Soc. 2016 (in print)
module lorenz96
    use params
    use utils, only: sum_2d
    use rp_emulator

    implicit none

    private

    !===========================================================================
    ! Interfaces
    !===========================================================================

    public step

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

    !===========================================================================
    ! Model parameters
    !===========================================================================
    
    public f, h, c, b, e, d, g_Z
    real(dp), parameter :: f = 20._dp
    real(dp), parameter :: h = 1._dp
    real(dp), parameter :: c = 10._dp
    real(dp), parameter :: b = 10._dp
    real(dp), parameter :: e = 10._dp
    real(dp), parameter :: d = 10._dp
    real(dp), parameter :: g_Z = 1.0_dp

    contains
        !=======================================================================
        ! Full three-level model
        !=======================================================================

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Steps forward once using 4th order Runge-Kutta scheme (full
        !> three-level model).
        !> @param[in] prev the previous state
        !> @return step the next state
         pure function step(prev)
            real(dp), dimension(truth_dim), intent(in) :: prev
            real(dp), dimension(truth_dim) :: step, k1, k2, k3, k4

            ! 4th order Runge-Kutta
            k1 = ode(prev)
            k2 = ode(prev+0.5_dp*dt*k1)
            k3 = ode(prev+0.5_dp*dt*k2)
            k4 = ode(prev+dt*k3) 

            step = prev + (dt/6._dp)*(k1 + 2._dp*k2 + 2._dp*k3 + k4)
        end function
        
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> The full three-level system of ODEs.
        !> @param[in] state state vector of X, Y and Z variables
        !> @return ode the evaluated ODE
        pure function ode(state)
            real(dp), intent(in) :: state(truth_dim)
            real(dp) :: x(n_x)
            real(dp) :: y(n_x*n_y)
            real(dp) :: z(n_x*n_y*n_z)
            real(dp) :: ode(truth_dim)

            ! Break up state vector into components
            x = state(:n_x)
            y = state(n_x+1:n_x+n_x*n_y)
            z = state(n_x+n_x*n_y+1:)

            ! Find derivative of each component separately
            ode(:n_x) = dXdT(x, y)
            ode(n_x+1:n_x+n_x*n_y) = dYdT(x, y, z)
            ode(n_x+n_x*n_y+1:) = dZdT(y, z)
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> The ODE for the X variables.
        !> @param[in] x state vector of X variables
        !> @param[in] y state vector of Y variables
        !> @return dXdT the evaluated ODE
        {% for type in types %}
        pure function dXdT_{{ type.name }}(x, y) result(dXdT)
            {{ type.code }}, intent(in) :: x(n_x)
            {{ type.code }}, intent(in) :: y(n_x*n_y)
            {{ type.code }} :: dXdT(n_x)
            {{ type.code }} :: sum_y(n_x)
            {{ type.code }} :: h_, c_, b_, f_

            h_ = h
            c_ = c
            b_ = b
            f_ = f

            ! Sum all y's for each x, making an n_x length vector of y sums
            sum_y = sum_2d(reshape(y, (/n_y,n_x/)))

            dXdT = cshift(x,-1)*(cshift(x,1)-cshift(x,-2)) - x + f_ & 
                & - (h_*c_/b_)*sum_y
        end function
        {% endfor %}

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> The ODE for the Y variables.
        !> @param[in] x state vector of X variables
        !> @param[in] y state vector of Y variables
        !> @param[in] z state vector of Z variables
        !> @return dYdT the evaluated ODE
        pure function dYdT(x, y, z)
            real(dp), intent(in) :: x(n_x)
            real(dp), intent(in) :: y(n_x*n_y)
            real(dp), intent(in) :: z(n_x*n_y*n_z)
            real(dp) :: dYdT(n_x*n_y)
            real(dp) :: x_rpt(n_x*n_y)
            real(dp) :: sum_z(n_x*n_y)
            integer :: k

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            ! Sum all z's for each y, making an n_x*n_y length vector of z sums
            sum_z = sum(reshape(z, (/n_z,n_x*n_y/)), dim=1)

            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - c*y
            dYdT = dYdT + (h*c/b)*x_rpt
            dYdT = dYdT - (h*e/d)*sum_z
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> The ODE for the Z variables.
        !> @param[in] y state vector of Y variables
        !> @param[in] z state vector of Z variables
        !> @return dZdT the evaluated ODE
        pure function dZdT(y, z)
            real(dp), intent(in) :: y(n_x*n_y)
            real(dp), intent(in) :: z(n_x*n_y*n_z)
            real(dp) :: dZdT(n_x*n_y*n_z)
            real(dp) :: y_rpt(n_x*n_y*n_z)
            integer :: k

            ! Repeat elements of y n_z times
            y_rpt = (/ (y(1+(k-1)/n_z), k = 1, n_x*n_y*n_z) /)

            dZdT = e*d*cshift(z,-1)*(cshift(z,1)-cshift(z,-2)) - g_Z*e*Z
            dZdT = dZdT + (h*e/d)*y_rpt
        end function
        
        !=======================================================================
        ! Three-level model with parametrised Z dynamics
        !=======================================================================     

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Steps forward once using 4th order Runge-Kutta scheme (two-level
        !> model with paramtrized Z level).
        !> @param[in] prev the previous state
        !> @param[in] stoch a vector of AR(1) random numbers used for the
        !> additive stochastic parametrization scheme.
        !> @return step the next state
        {% for type in types %}
        pure function step_param_z_{{ type.name }}(prev, stoch) result(step)
            {{ type.code }}, dimension(n_x+n_x*n_y), intent(in) :: prev
            {{ type.code }}, dimension(n_x*n_y), intent(in) :: stoch
            {{ type.code }}, dimension(n_x+n_x*n_y) :: step, k1, k2, k3, k4  
            {{ type.code }} :: half, six, two, dt_rpe

            half = 0.5_dp
            six = 6.0_dp
            two = 2.0_dp
            dt_rpe = dt
            
            ! 4th order Runge-Kutta
            k1 = ode_param_z(prev, stoch)
            k2 = ode_param_z(prev+half*dt_rpe*k1, stoch)
            k3 = ode_param_z(prev+half*dt_rpe*k2, stoch)
            k4 = ode_param_z(prev+dt_rpe*k3, stoch)

            step = prev + (dt_rpe/six)*(k1 + two*k2 + two*k3 + k4)
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> The two-level system of ODEs for the Lorenz '96 system, with
        !> parametrized Z.
        !> @param[in] state state vector of X and Y variables
        !> @param[in] stoch a vector of AR(1) random numbers used for the
        !> additive stochastic parametrization scheme.
        !> @return ode the evaluated ODE
        pure function ode_param_z_{{ type.name }}(state, stoch) result(ode)
            {{ type.code }}, intent(in) :: state(state_dim)
            {{ type.code }}, intent(in) :: stoch(n_x*n_y)
            {{ type.code }} :: x(n_x)
            {{ type.code }} :: y(n_x*n_y)
            {{ type.code }} :: ode(state_dim)

            ! Break up state vector into components
            x = state(:n_x)
            y = state(n_x+1:)

            ! Find derivative of each component separately
            ode(:n_x) = dXdT(x, y)
            ode(n_x+1:) = dYdT_param_z(x, y, stoch)
        end function
        
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> The ODE for the Y variables, with a parametrized Z-level.
        !> @param[in] x state vector of X variables
        !> @param[in] y state vector of Y variables
        !> @param[in] stoch a vector of AR(1) random numbers used for the
        !> additive stochastic parametrization scheme.
        !> @return dYdT the evaluated ODE
        pure function dYdT_param_z_{{ type.name }}(x, y, stoch) result(dYdT)
            {{ type.code }}, intent(in) :: x(n_x)
            {{ type.code }}, dimension(n_x*n_y), intent(in) :: y, stoch
            {{ type.code }}, dimension(n_x*n_y) :: dYdT, x_rpt, tend_z
            {{ type.code }}, dimension(5) :: coeffs
            {{ type.code }} :: h_, c_, b_
            integer :: k

            ! Set coefficients for deterministic part of paramtrization scheme
            ! (this is a 4th order polynomial)
            coeffs(1) = 0.001892_dp
            coeffs(2) = -0.066811_dp
            coeffs(3) = 0.131826_dp
            coeffs(4) = 0.242742_dp
            coeffs(5) = 0.039970_dp

            h_ = h
            c_ = c
            b_ = b

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            ! Compute the Z tendency from the chosen parametrization scheme
            tend_z = (coeffs(1)*y**4) + (coeffs(2)*y**3) + &
                & (coeffs(3)*y**2) + (coeffs(4)*y) + coeffs(5)
                
            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - c_*y
            dYdT = dYdT + (h_*c_/b_)*x_rpt        

            ! Apply additive stochastic term
            dYdT = dYdT - (tend_z + stoch)
        end function
        {% endfor %}
end module
