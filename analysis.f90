module analysis
    use params
    use utils

    implicit none

    public :: observe
    interface observe
        module procedure observe
        module procedure observe_rpe
    end interface observe

    contains
        pure function observe(state)
            real(dp), dimension(:, :), intent(in) :: state
            real(dp), dimension(obs_dim, size(state, 2)) :: observe

            observe = state(n_x+1:n_x+n_x*n_y,:)
        end function observe
        
        pure function observe_rpe(state)
            type(rpe_var), dimension(:, :), intent(in) :: state
            type(rpe_var), dimension(obs_dim, size(state, 2)) :: observe_rpe

            observe_rpe = observe(state%val)
        end function observe_rpe

        function enkf_assimilate(ensemble, obs_vec, obs_covar) result(analy)
            DOUBLE_OR_RPE, dimension(state_dim, n_ens), intent(in) :: ensemble
            real(dp), dimension(obs_dim), intent(in) :: obs_vec
            real(dp), dimension(obs_dim, obs_dim), intent(in) :: obs_covar
            DOUBLE_OR_RPE, dimension(state_dim, n_ens) :: analy, A
            DOUBLE_OR_RPE, dimension(state_dim) :: ens_mean
            real(dp), dimension(obs_dim, n_ens) ::  obs_table
            DOUBLE_OR_RPE, dimension(state_dim, obs_dim) :: gain
            DOUBLE_OR_RPE, dimension(state_dim, obs_dim) :: ens_cov_h_t
            real(dp) :: y_std
            integer :: i, j

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
                A(:, i) = ensemble(:, i) - ens_mean
            end do

            ! Ensemble covariance times transpose of observation matrix
            ens_cov_h_t = (rho/real(n_ens-1, dp)) * matmul(A, transpose(observe(A)))

            ! Kalman gain
            gain = matmul(ens_cov_h_t, inv(observe(ens_cov_h_t) + obs_covar))

            analy = ensemble + matmul(gain, obs_table - observe(ensemble))
        end function enkf_assimilate
end module analysis
