!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains various utility functions, e.g. for calculating statistics, sums.
!> Jinja2 templates are used to generate overloaded functions of different
!> precisions.
module utils
    use params
    use rp_emulator

    implicit none

    private

    !===========================================================================
    ! Overloaded functions
    !=========================================================================== 

    public std, identity

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

    public :: rmse_ens_mean
    interface rmse_ens_mean
        {% for type in types %}
        module procedure rmse_ens_mean_{{ type.name }}
        {% endfor %}
    end interface rmse_ens_mean

    public :: mean_ens_rmse
    interface mean_ens_rmse
        {% for type in types %}
        module procedure mean_ens_rmse_{{ type.name }}
        {% endfor %}
    end interface mean_ens_rmse

    public :: real
    interface real
        module procedure rpe_to_real
    end interface

    !===========================================================================
    ! Function definitions
    !=========================================================================== 

    contains
        {% for type in types %}
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Generates a random number drawn for the specified normal
        !> distribution.
        !> @param mean the mean of the distribution to draw from
        !> @param stdev the standard deviation of the distribution to draw from
        !> @return randn the generated random number
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

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Generates a random number drawn from an AR(1) process. This has been
        !> tuned for parametrizing the Z variables in the two-layer Lorenz '96
        !> model. See Wilkes, "Statistical Methods in the Atmospheric
        !> Sciences", 3rd ed. p. 410.
        !> @param last the last value of the process
        !> @return e the next value of the process
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

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Sum of 2D array (returns 1D array, with sum of elements along 1st
        !> dimension).
        !> @param array the input array
        !> @return sum_2d the sum
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

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> RMSE between ensemble mean and truth. Only computes RMSE across
        !> X-variables.
        !> @param ensemble the ensemble
        !> @param truth the truth state vector
        !> @return rmse_ens_mean the RMSE of the ensemble mean
        function rmse_ens_mean_{{ type.name }}(ensemble, truth) result(rmse_ens_mean)
            {{ type.code }}, dimension(state_dim, n_ens) :: ensemble
            real(dp), dimension(truth_dim) :: truth
            real(dp), dimension(n_x) :: ens_mean
            real(dp) :: rmse_ens_mean
            integer :: i

            ens_mean = (/ (sum(ensemble(i, :))/real(n_ens) , i = 1, n_x) /)

            rmse_ens_mean = sqrt(sum((ens_mean - truth(:n_x))**2)/real(n_x))
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Average RMSE of each ensemble member. Only computed over X-variables.
        !> @param ensemble the ensemble
        !> @param truth the truth state vector
        !> @return mean_ens_rmse the mean of the RMSE of each ensemble member
        function mean_ens_rmse_{{ type.name }}(ensemble, truth) result(mean_ens_rmse)
            {{ type.code }}, dimension(state_dim, n_ens) :: ensemble
            real(dp) :: truth(truth_dim)
            real(dp) :: ens_mse(n_ens)
            real(dp) :: mean_ens_rmse
            integer :: i

            ens_mse = (/ ( &
                sum((ensemble(:n_x,i) - truth(:n_x))**2)/real(n_x), i = 1, n_ens &
            & ) /)

            mean_ens_rmse = sum(sqrt(ens_mse))/real(n_ens,dp)
        end function
        {% endfor %}

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Calculates (biased) standard deviation of list of input variables.
        !> @param vars the input data
        !> @return std the standard deviation
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

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Generates an identity matrix of the given size.
        !> @param dimen the dimension of the required identity matrix
        !> @return identity the identity matrix (dimen x dimen)
        function identity(dimen)
            integer, intent(in) :: dimen
            real(dp) :: identity(dimen, dimen)
            integer :: i, j

            forall(i = 1:dimen, j = 1:dimen) identity(i, j) = real((i/j)*(j/i), dp)
        end function
 
        !===========================================================================
        ! Reduced precision utils
        !=========================================================================== 
        
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> A reduced precision matrix multiplier.
        !> @param mat1 the first matrix
        !> @param mat2 the second matrix
        !> @return matmul_rpe the multiplied matrix
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

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> A reduced precision matrix-vector multiplier.
        !> @param mat1 the matrix
        !> @param mat2 the vector
        !> @return matmul_rpe the vector pre-multiplied by the matrix
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
        
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Casts an rpe variable to a real value.
        !> @param rpe_input the input rpe variable
        !> @return rpe_to_real the double precision real value
        elemental real(dp) pure function rpe_to_real(rpe_input)
            type(rpe_var), intent(in) :: rpe_input

            rpe_to_real = rpe_input%val
        end function rpe_to_real
end module
