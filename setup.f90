!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains various functions and subroutines for setting up the assimilation.
module setup
    use params
    use lorenz96
    use utils, only: randn
    use rp_emulator

    implicit none

    private
    
    public spin_up, gen_ensemble, time_seed

    contains
        !> @brief
        !> Spins up the truth run.
        !> @return truth the spun-up state vector used for generating the truth
        !> run
        function spin_up() result(truth)
            real(dp) :: truth(truth_dim)
            integer :: i

            ! Some random initial conditions (doesn't really matter)
            truth(:n_x) = (/ (8, i = 1, n_x) /)
            truth(n_x+1:n_x+n_x*n_y) = (/ (randn(0._dp, 0.5_dp), i = 1, n_x*n_y) /)
            truth(n_x+n_x*n_y+1:) = (/ (randn(0._dp, 0.5_dp), i = 1, n_x*n_y*n_z) /)
            truth(4) = 8.008_dp

            ! Spin up
            do i = 1, 5000
                truth = step(truth)
            end do
        end function

        !> @brief
        !> Generate the first background ensemble by sampling from the truth
        !> run.
        !> @param truth the truth run to sample from
        !> @return ensemble the generated ensemble
        function gen_ensemble(truth) result(ensemble)
            real(dp), intent(in) :: truth(state_dim, n_steps)
            PRECISION :: ensemble(state_dim, n_ens)
            real(dp) :: rand
            integer :: i, j

            do i = 1, n_ens
                call random_number(rand)
                j = ceiling(rand * n_steps)
                ensemble(:, i) = truth(:, j)
            end do
        end function

        !> @brief
        !> Seeds RNG from system clock.
        subroutine time_seed()
          integer :: i, n, clock
          integer, allocatable :: seed(:)
        
          call random_seed(size = n)
          allocate(seed(n))
        
          call system_clock(count=clock)
        
          seed = clock + 37 * (/ (i - 1, i = 1, n) /)
          call random_seed(put = seed)
        
          deallocate(seed)
        end subroutine
end module
