module utils
    use params, only: dp, n_x, n_y

    implicit none

    contains
        ! Generates a random number drawn for the specified normal distribution
        function randn(mean, stdev)
            real(dp), intent(in) :: mean, stdev
            real(dp) :: randn, rand(2), u, v

            call random_number(rand)

            ! Box-Muller method
            u = (-2.0_dp * log(rand(1))) ** 0.5_dp
            v =   2.0_dp * 6.28318530718_dp * rand(2)
            randn = mean + stdev * u * sin(v)
        end function randn

        ! Seeds from system clock
        subroutine time_seed()
          integer :: i, n, clock
          integer, dimension(:), allocatable :: seed
        
          call random_seed(size = n)
          allocate(seed(n))
        
          call system_clock(count=clock)
        
          seed = clock + 37 * (/ (i - 1, i = 1, n) /)
          call random_seed(put = seed)
        
          deallocate(seed)
        end subroutine

        ! Matrix inverter based on LU decomposition
        ! Depends on LAPACK
        function inv(m) result(m_inv)
            real(dp), dimension(:,:), intent(in) :: m
            real(dp), dimension(size(m, 1), size(m, 2)) :: m_inv
            real(dp), dimension(size(m, 1)) :: work
            integer, dimension(size(m, 1)) :: ipiv
            integer :: n, info

            external rgetrf
            external rgetri

            m_inv = m
            n = size(m, 1)

            call rgetrf(n, m_inv, n, ipiv, info)

            if (info > 0) then
                stop 'Singular matrix!'
            else if (info < 0) then
                stop 'Illegal argument to dgetrf'
            end if

            call rgetri(n, m_inv, n, ipiv, work, n, info)

            if (info /= 0) then
                stop 'Matrix inversion failed'
            end if
        end function inv
        
        ! Calculate (biased) standard deviation
        function std(vars)
            real(dp), intent(in) :: vars(:)
            real(dp) :: std
            real(dp) :: mean
            integer :: n, j
            n = size(vars)
            
            mean = sum(vars) / real(n, dp)
            std = sum((/ ((vars(j) - mean)**2, j = 1, n) /))
            std = sqrt(std / real(n, dp))
        end function
        
        ! Zero mean AR(1) process
        function additive_noise(last) result(e)
            real(dp), dimension(n_x*n_y), intent(in) :: last
            real(dp), dimension(n_x*n_y) :: e
            real(dp) :: phi = 0.997_dp
            real(dp) :: sigma_e = 0.126_dp
            real(dp), dimension(n_x*n_y) :: z
            integer :: i
            
            do i = 1, n_x*n_y
                z(i) = randn(0.0_dp, sigma_e)
            end do
            
            e = phi * last + sqrt(1-phi**2) * z
        end function additive_noise
end module utils
