module rp_utils
    use rp_emulator

    implicit none
    
    public :: inv
    interface inv
        module procedure inv_rpe
    end interface inv 
    
    contains
        ! Matrix inverter based on LU decomposition
        ! Depends on LAPACK
        function inv_rpe(m) result(m_inv)
            type(rpe_var), dimension(:,:), intent(in) :: m
            type(rpe_var), dimension(size(m, 1), size(m, 2)) :: m_inv
            type(rpe_var), dimension(size(m, 1)) :: work
            integer, dimension(size(m, 1)) :: ipiv
            integer :: n, info

            external rgetrf
            external rgetri

            m_inv = m
            n = size(m, 1)

            call rgetrf(n, n, m_inv, n, ipiv, info)

            if (info > 0) then
                stop 'Singular matrix!'
            else if (info < 0) then
                stop 'Illegal argument to dgetrf'
            end if

            call rgetri(n, m_inv, n, ipiv, work, n, info)

            if (info /= 0) then
                stop 'Matrix inversion failed'
            end if
        end function inv_rpe
end module rp_utils
