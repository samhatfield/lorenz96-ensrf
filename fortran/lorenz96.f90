module lorenz96
    use params

    implicit none

    ! Model parameters
    real(dp), parameter :: f = 20
    real(dp) :: h = 1
    real(dp) :: c = 4
    real(dp) :: b = 10

    contains
        ! Step forward once
        pure function step(prev_state)
            real(dp), dimension(state_dim), intent(in) :: prev_state
            real(dp), dimension(n_x) :: x, k1, k2, k3, k4
            real(dp), dimension(n_x*n_y) :: y, l1, l2, l3, l4
            real(dp), dimension(state_dim) :: step

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
            real(dp), dimension(n_x), intent(in) :: x
            real(dp), dimension(n_x*n_y), intent(in) :: y
            real(dp), dimension(n_x) :: dXdT
            real(dp), dimension(n_x) :: sum_y

            sum_y = sum(reshape(y, (/n_y,n_x/)), dim=1)

            dXdT = cshift(x,-1)*(cshift(x,1)-cshift(x,-2)) - x + f
            dXdT = dXdT - (h*c/b)*sum_y
        end function dXdT

        ! Y ODE
        pure function dYdT(x, y)
            real(dp), dimension(n_x), intent(in) :: x
            real(dp), dimension(n_x*n_y), intent(in) :: y
            real(dp), dimension(n_x*n_y) :: dYdT
            real(dp), dimension(n_x*n_y) :: x_rpt
            integer :: k

            ! Repeat elements of x n_y times
            x_rpt = (/ (x(1+(k-1)/n_y), k = 1, n_x*n_y) /)

            dYdT = c*b*cshift(y,1)*(cshift(y,-1)-cshift(y,2)) - c*y
            dYdT = dYdT + (h*c/b)*x_rpt
        end function dYdT
end module lorenz96
