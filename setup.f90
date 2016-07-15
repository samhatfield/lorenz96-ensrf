!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains various functions and subroutines for setting up the assimilation.
module setup
    use params
    use lorenz96
    use utils, only: randn

    implicit none

    private
    
    public spin_up, gen_ensemble, write_params, time_seed

    contains
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
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

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Generate the first background ensemble by sampling from the truth
        !> run.
        !> @param truth the truth run to sample from
        !> @return ensemble the generated ensemble
        function gen_ensemble(truth) result(ensemble)
            real(dp), intent(in) :: truth(truth_dim, n_steps)
            PRECISION :: ensemble(state_dim, n_ens)
            real(dp) :: rand
            integer :: i, j

            do i = 1, n_ens
                call random_number(rand)
                j = ceiling(rand * n_steps)
                ensemble(:, i) = truth(:n_x+n_x*n_y, j)
            end do
        end function

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
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

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Writes the header of the output YAML file. This is basically just a
        !> list of key-value pairs giving the parameters for this run.
        subroutine write_params()
            integer :: file_1 = 20

            open(unit=file_1, file="results.yml", action="write", status="replace")
            1 format(a4, f5.3)
            2 format(a4, i3)
            3 format(a6, i4)
            4 format(a11, i3)
            5 format(a8, f5.2)
            6 format(a4, f7.2)
            7 format(a3, f6.2)
            8 format(a8, i5)
            9 format(a6, f7.4)

            write (file_1, 1) 'dt: ', dt 
            write (file_1, 6) 'fin: ', fin
            write (file_1, 2) 'n_x: ', n_x
            write (file_1, 2) 'n_y: ', n_y
            write (file_1, 2) 'n_z: ', n_z
            write (file_1, 3) 'n_ens: ', n_ens
            write (file_1, 7) 'h: ', h
            write (file_1, 7) 'c: ', c
            write (file_1, 7) 'b: ', b
            write (file_1, 7) 'e: ', e
            write (file_1, 7) 'd: ', d
            write (file_1, 6) 'g_Z: ', g_Z 
            write (file_1, 4) 'assim_freq: ', assim_freq
            write (file_1, 9), 'y_var: ', y_var
            write (file_1, 8) 'obs_dim: ', obs_dim
            write (file_1, 3) 'sbits: ', sbits
            write (file_1, '(a)') '# Columns (all quantities calculated from &
            & X variables only:'
            write (file_1, '(a)') '# Avg norm, Std norm, Truth'
            write (file_1, '(a5)') '---'

            close(file_1)
        end subroutine write_params
end module
