module utils
    use params, only: dp, n_x, n_y
    use rp_emulator

    implicit none

    !===========================================================================
    ! Overloaded functions
    !=========================================================================== 

    public :: inv
    interface inv
        module procedure inv
        module procedure inv_rpe
    end interface inv
    
    public :: randn
    interface randn
        module procedure randn
        module procedure randn_rpe
    end interface randn
    
    public :: ar_1
    interface ar_1
        module procedure ar_1
        module procedure ar_1_rpe
    end interface ar_1

    public :: matmul
    interface matmul
        module procedure matmul_rpe
    end interface matmul

    public :: sum_1d
    interface sum_1d
        module procedure sum_1d
        module procedure sum_1d_rpe
    end interface sum_1d

    public :: real
    interface real
        module procedure rpe_to_real
    end interface real

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
        function ar_1(last) result(e)
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
        end function ar_1

        pure function sum_1d(array)
        	real(dp), intent(in) :: array(:)
        	real(dp) :: sum_1d
        	integer :: i, n
        	
        	n = size(array)
        	
        	sum_1d = 0.0d0
        	
        	do i = 1, n
        		sum_1d = sum_1d + array(i)
        	end do
        end function sum_1d
 
        !===========================================================================
        ! Reduced precision utils
        !=========================================================================== 
        
        ! Generates a random number drawn for the specified normal distribution
        function randn_rpe(mean, stdev)
            type(rpe_var), intent(in) :: mean, stdev
            type(rpe_var) :: randn_rpe
            
            randn_rpe = randn(mean%val, stdev%val)
        end function randn_rpe
        
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
        
        ! Zero mean AR(1) process
        function ar_1_rpe(last) result(e)
            type(rpe_var), dimension(n_x*n_y), intent(in) :: last
            type(rpe_var), dimension(n_x*n_y) :: e
            real(dp) :: phi = 0.997_dp
            real(dp) :: sigma_e = 0.126_dp
            real(dp), dimension(n_x*n_y) :: z
            integer :: i
            
            do i = 1, n_x*n_y
                z(i) = randn(0.0_dp, sigma_e)
            end do
            
            e = phi * last + sqrt(1-phi**2) * z
        end function ar_1_rpe
        
        ! Reduced precision matrix multiplication
        function matmul_rpe(mat1, mat2)
            type(rpe_var), intent(in) :: mat1(:,:), mat2(:,:)
            type(rpe_var) :: matmul_rpe(size(mat1, 1), size(mat2, 2))
            integer :: i, j, k
            
            matmul_rpe(:,:) = 0.0d0
            
            do i = 1, size(mat2, 2)
                do j = 1, size(mat1, 1)
                    do k = 1, size(mat1, 2)
                        matmul_rpe(j, i) = matmul_rpe(j, i) + mat1(j, k)*mat2(k, i)
                    end do
                end do
            end do
        end function matmul_rpe
        
        ! Reduced precision sum of 2D array (returns 1D array, with sum of elements
        ! along 1st dimension)
        pure function sum_2d_rpe(array)
        	type(rpe_var), intent(in) :: array(:,:)
        	type(rpe_var) :: sum_2d_rpe(size(array, 1))
        	integer :: i, n
        	
        	n = size(array, 1)
        	
        	sum_2d_rpe(:) = 0.0d0
        	
        	do i = 1, n
        		sum_2d_rpe = sum_2d_rpe + array(i, :)
        	end do
        end function sum_2d_rpe
        
        pure function sum_1d_rpe(array)
        	type(rpe_var), intent(in) :: array(:)
        	type(rpe_var) :: sum_1d_rpe
        	integer :: i, n
        	
        	n = size(array)
        	
        	sum_1d_rpe = 0.0d0
        	
        	do i = 1, n
        		sum_1d_rpe = sum_1d_rpe + array(i)
        	end do
        end function sum_1d_rpe

        elemental real pure function rpe_to_real(rpe_input)
            type(rpe_var), intent(in) :: rpe_input

            rpe_to_real = rpe_input%val
        end function rpe_to_real
end module utils
