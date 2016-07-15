module utils
    use params
    use rp_emulator

    implicit none

    !===========================================================================
    ! Overloaded functions
    !=========================================================================== 

    public :: randn
    interface randn
        {% for type in types %}
        module procedure randn_{{ type.name }}
        {% endfor %}
    end interface
    
    public :: ar_1
    interface ar_1
        {% for type in types %}
        module procedure ar_1_{{ type.name }}
        {% endfor %}
    end interface

    public :: sum_2d
    interface sum_2d
        {% for type in types %}
        module procedure sum_2d_{{ type.name }}
        {% endfor %}
    end interface

    public :: matmul
    interface matmul
        module procedure matmul_rpe
        module procedure matmulvec_rpe
    end interface

    public :: sum_1d
    interface sum_1d
        {% for type in types %}
        module procedure sum_1d_{{ type.name }}
        {% endfor %}
    end interface

    public :: rmse_mean
    interface rmse_mean
        {% for type in types %}
        module procedure rmse_mean_{{ type.name }}
        {% endfor %}
    end interface rmse_mean

    public :: real
    interface real
        module procedure rpe_to_real
    end interface

    contains
        ! Generates a random number drawn for the specified normal distribution
        {% for type in types %}
        function randn_{{ type.name }}(mean, stdev) result(randn)
            {{ type.code }}, intent(in) :: mean, stdev
            {{ type.code }} :: u, v, randn
            real(dp) :: rand(2)

            call random_number(rand)

            ! Box-Muller method
            u = (-2.0_dp * log(rand(1))) ** 0.5_dp
            v =   2.0_dp * 6.28318530718_dp * rand(2)
            randn = mean + stdev * u * sin(v)
        end function

        ! Zero mean AR(1) process
        function ar_1_{{ type.name }}(last) result(e)
            {{ type.code }}, intent(in) :: last(n_x*n_y)
            {{ type.code }} :: e(n_x*n_y)
            {{ type.code }} :: phi, sigma_e, zero
            {{ type.code }} :: z(n_x*n_y)
            integer :: i

            phi = 0.997_dp
            sigma_e = 0.126_dp
            zero = 0.0_dp
            
            do i = 1, n_x*n_y
                z(i) = randn(zero, sigma_e)
            end do
            
            e = phi * last + sqrt(1-phi**2) * z
        end function

        ! Reduced precision sum of 2D array (returns 1D array, with sum of elements
        ! along 1st dimension)
        pure function sum_2d_{{ type.name }}(array) result(sum_2d)
            {{ type.code }}, intent(in) :: array(:,:)
        	{{ type.code }} :: sum_2d(size(array, 1))
        	integer :: i, n
        	
        	n = size(array, 1)
        	
        	sum_2d(:) = 0.0_dp
        	
        	do i = 1, n
        		sum_2d = sum_2d + array(i, :)
        	end do
        end function

        pure function sum_1d_{{ type.name }}(array) result(sum_1d)
            {{ type.code }}, intent(in) :: array(:)
            {{ type.code }} :: sum_1d
        	integer :: i, n
        	
        	n = size(array)
        	
        	sum_1d = 0.0_dp
        	
        	do i = 1, n
        		sum_1d = sum_1d + array(i)
        	end do
        end function

        function rmse_mean_{{ type.name }}(ensemble, truth) result(rmse_mean)
            {{ type.code }}, dimension(state_dim, n_ens) :: ensemble
            real(dp), dimension(truth_dim) :: truth
            real(dp), dimension(n_x) :: ens_mean
            real(dp) :: rmse_mean
            integer :: i

            rmse_mean = 0.0_dp

            ens_mean = (/ (sum_1d(ensemble(i, :))/real(n_ens) , i = 1, n_x) /)

            rmse_mean = sqrt(sum((ens_mean - truth(:n_x))**2)/real(n_x,dp))
        end function
        {% endfor %}

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

        ! Generate identity matrix
        function identity(dimen)
            integer, intent(in) :: dimen
            real(dp) :: identity(dimen, dimen)
            integer :: i, j

            forall(i = 1:dimen, j = 1:dimen) identity(i, j) = real((i/j)*(j/i), dp)
        end function
 
        !===========================================================================
        ! Reduced precision utils
        !=========================================================================== 
        
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

        ! Reduced precision matrix-vector multiplication
        function matmulvec_rpe(mat, vec)
            type(rpe_var), intent(in) :: mat(:,:), vec(:)
            type(rpe_var) :: matmulvec_rpe(size(mat, 1))
            integer :: i, j, k
            
            do j = 1, size(mat, 1)
                matmulvec_rpe(j) = 0.0_dp
                do k = 1, size(vec)
                    matmulvec_rpe(j) = matmulvec_rpe(j) + mat(j, k)*vec(k)
                end do
            end do
        end function matmulvec_rpe
        
        elemental real pure function rpe_to_real(rpe_input)
            type(rpe_var), intent(in) :: rpe_input

            rpe_to_real = rpe_input%val
        end function rpe_to_real
end module utils
