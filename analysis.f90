!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions for performing ensemble square-root filtering and
!> localisation of gain matrices, with the Lorenz '96/'95 system in mind.
!>
!> Bibliography
!> ============
!> Whitaker and Hamill, "Ensemble Data Assimilation without Perturbed
!> Observations", Mon. Weather. Rev. 130 (2002)
!> Gaspari and Cohn, "Construction of correlation functions in two and three
!> dimensions", Q. J. R. Meterol. Soc 125 (1999)
module analysis
    use params
    use utils, only: matmul
    use observation, only: observe
    use rp_emulator

    implicit none

    private
    public ensrf_assimilate

    contains
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Performs sequential ensemble square-root filtering data assimilation
        !> with the given background ensemble, observation vector and
        !> observation covariance matrix.
        !> @param[in] background the background ensemble
        !> @param[in] obs the observation vector
        !> @param[in] R the observation covariance matrix (assumed diagonal)
        !> @return analysis the analysis ensemble
        function ensrf_assimilate(background, obs, R) result(analysis)
            PRECISION, intent(in) :: background(state_dim, n_ens)
            real(dp), intent(in) :: obs(obs_dim)
            real(dp), intent(in) :: R(obs_dim, obs_dim)
            PRECISION, dimension(state_dim, n_ens) :: analysis, X_f
            PRECISION, dimension(state_dim) :: P_f_H_T, gain, ens_mean
            PRECISION :: alpha, HP_f_H_T, one, rho_
            integer :: i, j

            one = 1.0_dp
            rho_ = rho

            ! Mean ensemble vector
            ens_mean = (/ (sum(background(j, :))/real(n_ens) , j = 1, state_dim) /)

            ! Form the background ensemble perturbation matrix (with covariance
            ! inflation)
            do i = 1, n_ens
                X_f(:, i) = rho_ * (background(:, i) - ens_mean)
            end do

            ! Sequentially process observations
            do i = 1, obs_dim
                ! Ensemble covariance times transpose of observation matrix
                P_f_H_T = matmul(X_f, observe(X_f, i)) / (n_ens-1)

                HP_f_H_T = observe(P_f_H_T, i)
    
                ! Kalman gain
                gain = P_f_H_T / (HP_f_H_T + R(i, i))

                ! Localize X variables
                gain(:n_x) = localize_x(gain(:n_x), y_skip*i-(y_skip-1))

                ! Localize Y variables
                gain(n_x+1:n_x+n_x*n_y) = localize_y(gain(n_x+1:n_x+n_x*n_y), y_skip*i-(y_skip-1))

                ! Update ensemble mean
                ens_mean = ens_mean + gain * (obs(i) - observe(ens_mean, i))

                ! Update perturbations
                alpha = one/(one+sqrt(R(i,i)/(HP_f_H_T+R(i,i))))
                do j = 1, n_ens
                    X_f(:, j) = X_f(:, j) - alpha*gain * observe(X_f(:, j), i)
                end do
            end do

            ! Form final ensemble
            do i = 1, n_ens
                analysis(:, i) = ens_mean + X_f(:, i)
            end do
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Localises given Y-variable Kalman gain vector based on the given
        !> currently assimilated Y variable.
        !> @param[in] gain the Kalman gain
        !> @param[in] y the Y variable currently being assimilated
        !> @return local_gain the localised gain vector
        function localize_y(gain, y) result(local_gain)
            PRECISION, intent(in) :: gain(n_x*n_y)
            integer, intent(in) :: y
            PRECISION :: local_gain(n_x*n_y)
            PRECISION, save :: mask(n_x*n_y)
            PRECISION :: c
            logical, save :: first_call = .True.
            integer :: i, j

            c = loc_y

            ! The localisation mask is only computed once, because it never
            ! changes
            if (first_call) then
                do i = 0, n_x*n_y-1
                    ! Account for cyclicity of Y variables
                    if (i >= n_x*n_y/2) then
                        j = n_x*n_y - i
                    else
                        j = i
                    end if
    
                    ! Generate mask based on Gaspari-Cohn function
                    mask(i+1) = gaspari_cohn(j, c)
                end do

                first_call = .False.
            end if

            ! Schur product
            local_gain = gain * cshift(mask, 1-y)
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Localises given X-variable Kalman gain vector based on the given
        !> currently assimilated Y variable.
        !> @param[in] gain the Kalman gain
        !> @param[in] y the Y variable currently being assimilated
        !> @return local_gain the localised gain vector
        function localize_x(gain, y) result(local_gain)
            PRECISION, intent(in) :: gain(n_x)
            integer, intent(in) :: y
            PRECISION :: local_gain(n_x)
            PRECISION, save :: mask(n_x)
            PRECISION :: c
            logical, save :: first_call = .True.
            integer :: i, j
            integer :: dist_y

            ! Determine the X-variable linked to the currently assimilated
            ! Y variable. The distance to this X variable is input to the
            ! localisation function
            dist_y = int((y-1)/n_y)+1

            c = loc_x

            ! The localisation mask is only computed once, because it never
            ! changes
            if (first_call) then
                do i = 0, n_x-1
                    ! Account for cyclicity
                    if (i >= n_x/2) then
                        j = n_x - i
                    else
                        j = i
                    end if

                    ! Generate mask based on Gaspari-Cohn function
                    mask(i+1) = gaspari_cohn(j, c)
                end do

                first_call = .False.
            end if

            ! Schur product
            local_gain = gain * cshift(mask, 1-dist_y)
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Gaspari-Cohn localisation function (see bibliography in file
        !> header).
        !> @param[in] z_ the distance between the two grid-points/variables
        !> @param[in] c the localisation length scale
        !> @return C_0 the value of the function
        pure function gaspari_cohn(z_, c) result(C_0)
            integer, intent(in) :: z_
            PRECISION, intent(in) :: c
            PRECISION :: z
            PRECISION :: C_0

            z = z_

            if (0 <= z .and. z <= c) then
                C_0 = -0.25_dp*(z/c)**5 + 0.5_dp*(z/c)**4 + 0.625_dp*(z/c)**3 - &
                    & (5.0_dp/3.0_dp)*(z/c)**2 + 1.0_dp
            else if (c < z .and. z <= 2*c) then
                C_0 = (1.0_dp/12.0_dp)*(z/c)**5 - 0.5_dp*(z/c)**4 + 0.625_dp*(z/c)**3 + &
                    & (5.0_dp/3.0_dp)*(z/c)**2 - 5.0_dp*(z/c) + 4.0_dp - (2.0_dp/3.0_dp)*(c/z)
            else
                C_0 = 0.0_dp
            end if
        end function
end module
