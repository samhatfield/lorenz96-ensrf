module analysis
    use params
    use utils

    implicit none

    contains
        pure function observe(state)
            real(dp), dimension(:, :), intent(in) :: state
            real(dp), dimension(obs_dim, size(state, 2)) :: observe

            observe = state
        end function observe

        function assimilate(ensemble, obs_vec, obs_covar, sig_obs) result(analy)
            real(dp), dimension(state_dim, n_ens), intent(in) :: ensemble
            real(dp), dimension(obs_dim), intent(in) :: obs_vec
            real(dp), dimension(obs_dim, obs_dim), intent(in) :: obs_covar
            real(dp), intent(in) :: sig_obs
            real(dp), dimension(state_dim, n_ens) :: analy, A
            real(dp), dimension(state_dim) :: ens_mean
            real(dp), dimension(obs_dim, n_ens) ::  obs_table
            real(dp), dimension(state_dim, obs_dim) :: gain
            real(dp), dimension(state_dim, obs_dim) :: ens_cov_h_t
            real(dp) :: rho = 1.4
            integer :: i, j

            ! Mean ensemble vector
            ens_mean = (/ (sum(ensemble(i, :))/real(n_ens) , i = 1, state_dim) /)

            ! Table of perturbed observations
            do i = 1, n_ens
                do j = 1, obs_dim
                    obs_table(j, i) = obs_vec(j) + randn(0.0d0, sig_obs)
                end do
            end do

            ! 'Anomaly table'
            do i = 1, n_ens
                A(:, i) = ensemble(:, i) - ens_mean
            end do

            ! Ensemble covariance times transpose of observation matrix
            ens_cov_h_t = (rho/(n_ens-1)) * matmul(A, transpose(observe(A)))

            ! Kalman gain
            gain = matmul(ens_cov_h_t, inv(observe(ens_cov_h_t) + obs_covar))

            analy = ensemble + matmul(gain, obs_table - observe(ensemble))
        end function assimilate
end module analysis
