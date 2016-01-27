module utils
    use params, only: dp

    implicit none

    contains
        ! Generates a random number drawn for the specified normal distribution
        function randn(mean, stdev)
            real(dp), intent(in) :: mean, stdev
            real(dp) :: randn, rand(2), u, v

            call random_number(rand)

            ! Box-Muller method
            u = (-2.0d0 * log(rand(1))) ** 0.5
            v =   2.0d0 * 6.28318530718 * rand(2)
            randn = mean + stdev * u * sin(v)
        end function randn

        ! Matrix inverter based on LU decomposition
        ! Depends on LAPACK
        function inv(m) result(m_inv)
            real(dp), dimension(:,:), intent(in) :: m
            real(dp), dimension(size(m, 1), size(m, 2)) :: m_inv
            real(dp), dimension(size(m, 1)) :: work
            integer, dimension(size(m, 1)) :: ipiv
            integer :: n, info

            external dgetrf
            external dgetri

            m_inv = m
            n = size(m, 1)

            call dgetrf(n, n, m_inv, n, ipiv, info)

            if (info > 0) then
                stop 'Singular matrix!'
            else if (info < 0) then
                stop 'Illegal argument to dgetrf'
            end if

            call dgetri(n, m_inv, n, ipiv, work, n, info)

            if (info /= 0) then
                stop 'Matrix inversion failed'
            end if
        end function inv
end module utils
