module lorenz96
    use params

    implicit none

    ! Model parameters
    real(dp), parameter :: f = 20
    real(dp), parameter :: h = 1
    real(dp), parameter :: c = 4
    real(dp), parameter :: b = 10
    real(dp), parameter :: e = 4
    real(dp), parameter :: d = 10
    real(dp), parameter :: g_X = 1
    real(dp), parameter :: g_Y = c
    real(dp), parameter :: g_Z = e

    contains
        ! Step forward once
        pure function step(prev_state)
            real(dp), dimension(state_dim), intent(in) :: prev_state
            real(dp), dimension(n_x) :: x, k1, k2, k3, k4
            real(dp), dimension(n_x*n_y) :: y, l1, l2, l3, l4
            real(dp), dimension(n_x*n_y*n_z) :: z, m1, m2, m3, m4
            real(dp), dimension(state_dim) :: step

            x = prev_state(:n_x)
            y = prev_state(n_x+1:n_x+n_x*n_y)
            z = prev_state(n_x+n_x*n_y+1:)
        
            ! 3 coupled ODEs 4th order Runge-Kutta
            k1 = dXdT(x, y)
            l1 = dYdT(x, y, z)
            m1 = dZdT(y, z)
            k2 = dXdT(x+0.5*dt*k1, y+0.5*dt*l1)
            l2 = dYdT(x+0.5*dt*k1, y+0.5*dt*l1, z+0.5*dt*m1)
            m2 = dZdT(y+0.5*dt*l1, z+0.5*dt*m1)
            k3 = dXdT(x+0.5*dt*k2, y+0.5*dt*l2)
            l3 = dYdT(x+0.5*dt*k2, y+0.5*dt*l2, z+0.5*dt*m2)
            m3 = dZdT(y+0.5*dt*l2, z+0.5*dt*m2)
            k4 = dXdT(x+dt*k3, y+dt*l3)
            l4 = dYdT(x+dt*k3, y+dt*l3, z+dt*m3)
            m4 = dZdT(y+dt*l3, z+dt*m3)
        
            x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
            y = y + (dt/6)*(l1 + 2*l2 + 2*l3 + l4)
            z = z + (dt/6)*(m1 + 2*m2 + 2*m3 + m4)

            ! Concatenate X's and Y's into return value
            step(:n_x) = x
            step(n_x+1:n_x+n_x*n_y) = y
            step(n_x+n_x*n_y+1:) = z
        end function step

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
end module lorenz96
