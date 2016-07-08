module analysis
    use params
    use utils
    use observation
    use local

    implicit none

    contains
        function ensrf_assimilate(background, obs, R) result(analysis)
            PRECISION, dimension(state_dim, n_ens), intent(in) :: background
            real(dp), dimension(obs_dim), intent(in) :: obs
            real(dp), dimension(obs_dim, obs_dim), intent(in) :: R
            PRECISION, dimension(state_dim, n_ens) :: analysis, X_f
            PRECISION, dimension(state_dim) :: ens_mean
            PRECISION, dimension(state_dim) :: P_f_H_T, gain
            PRECISION :: alpha, HP_f_H_T, one, rho_
            integer :: i, j

            one = 1.0_dp
            rho_ = rho

            ! Mean ensemble vector
            ens_mean = (/ (sum_1d(background(j, :))/real(n_ens) , j = 1, state_dim) /)

            ! Form the background ensemble perturbation matrix (with covariance inflation)
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

                ! Localize Y variables
                gain(n_x+1:n_x+n_x*n_y) = localize(gain(n_x+1:n_x+n_x*n_y), i)

                ! Update ensemble mean
                ens_mean = ens_mean + gain * (obs(i) - observe(ens_mean, i))

                ! Update pertubations
                alpha = one/(one+sqrt(R(i,i)/(HP_f_H_T+R(i,i))))
                do j = 1, n_ens
                    X_f(:, j) = X_f(:, j) - alpha*gain * observe(X_f(:, j), i)
                end do
            end do

            ! Form final ensemble
            do i = 1, n_ens
                analysis(:, i) = ens_mean + X_f(:, i)
            end do
        end function ensrf_assimilate
end module
