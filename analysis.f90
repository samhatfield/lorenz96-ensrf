module analysis
    use params
    use utils
    use observation

    implicit none

    contains
        function enkf_assimilate(ensemble, obs_vec, obs_covar) result(analy)
            DOUBLE_OR_RPE, dimension(state_dim, n_ens), intent(in) :: ensemble
            real(dp), dimension(obs_dim), intent(in) :: obs_vec
            real(dp), dimension(obs_dim, obs_dim), intent(in) :: obs_covar
            DOUBLE_OR_RPE, dimension(state_dim, n_ens) :: analy, X_f
            DOUBLE_OR_RPE, dimension(state_dim) :: ens_mean
            real(dp), dimension(obs_dim, n_ens) ::  obs_table
            DOUBLE_OR_RPE, dimension(state_dim, obs_dim) :: gain
            DOUBLE_OR_RPE, dimension(state_dim, obs_dim) :: ens_cov_h_t
            DOUBLE_OR_RPE :: one
            real(dp) :: y_std
            integer :: i, j

            one = 1.0_dp

            ! Get y standard deviation (assume all y variances in y covariance matrix are the same)
            y_std = sqrt(obs_covar(1, 1))

            ! Mean ensemble vector
            ens_mean = (/ (sum_1d(ensemble(i, :))/real(n_ens) , i = 1, state_dim) /)

            ! Table of perturbed observations
            do i = 1, n_ens
                do j = 1, obs_dim
                    obs_table(j, i) = obs_vec(j) + randn(0.0_dp, y_std)
                end do
            end do

            ! 'Anomaly table'
            do i = 1, n_ens
                X_f(:, i) = ensemble(:, i) - ens_mean
            end do

            ! Ensemble covariance times transpose of observation matrix
            ens_cov_h_t = (one/n_ens-1, dp) * matmul(X_f, transpose(observe(X_f)))

            ! Kalman gain
            gain = matmul(ens_cov_h_t, inv(observe(ens_cov_h_t) + obs_covar))

            analy = ensemble + matmul(gain, obs_table - observe(ensemble))
        end function enkf_assimilate

        function ensrf_assimilate(ensemble, obs_vec, R) result(analy)
            DOUBLE_OR_RPE, dimension(state_dim, n_ens), intent(in) :: ensemble
            real(dp), dimension(obs_dim), intent(in) :: obs_vec
            real(dp), dimension(obs_dim, obs_dim), intent(in) :: R
            DOUBLE_OR_RPE, dimension(state_dim, n_ens) :: analy, X_f
            DOUBLE_OR_RPE, dimension(state_dim) :: ens_mean
            DOUBLE_OR_RPE, dimension(state_dim) :: P_f_H_T, gain
            DOUBLE_OR_RPE :: alpha, HP_f_H_T, one, rho
            integer :: i, j

            one = 1.0_dp
            rho = 1.02_dp

            ! Mean ensemble vector
            ens_mean = (/ (sum_1d(ensemble(j, :))/real(n_ens) , j = 1, state_dim) /)

            ! Form the ensemble perturbation matrix
            do i = 1, n_ens
                X_f(:, i) = rho * (ensemble(:, i) - ens_mean)
            end do

            ! Sequentially process observations
            do i = 1, obs_dim
                ! Ensemble covariance times transpose of observation matrix
                P_f_H_T = (one/(n_ens-1))*matmul(X_f, observe(X_f, i))

                HP_f_H_T = observe_1d_row(P_f_H_T, i)
    
                ! Kalman gain
                gain = P_f_H_T / (HP_f_H_T + R(i, i))
    
                ! Update ensemble mean
                ens_mean = ens_mean + gain * (obs_vec(i) - observe_1d_row(ens_mean, i))

                ! Update pertubations
                alpha = one/(one+sqrt(R(i,i)/(HP_f_H_T+R(i,i))))
                do j = 1, n_ens
                    X_f(:, j) = X_f(:, j) - alpha*gain * observe_1d_row(X_f(:, j), i)
                end do
            end do

            ! Form final ensemble
            do i = 1, n_ens
                analy(:, i) = ens_mean + X_f(:, i)
            end do
        end function ensrf_assimilate
end module analysis
