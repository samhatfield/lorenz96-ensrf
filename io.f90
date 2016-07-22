!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions necessary for writing to an output NetCDF file.
module io
    use rp_emulator
    use params
    use netcdf
    use utils, only: real, rmse_mean

    implicit none

    private
    public setup_output, output, open_file, close_file

    integer :: ncid, timedim, timevar, xdim, xvar, truthx, ensdim, ensvar, ensx, rmsemeanx
    logical :: reduced_

    contains
        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Sets up the output NetCDF file, including writing simulation
        !> parameters.
        !> NOTE: uses a temporary file 'out.txt' to get the git revision.
        !> This file will be overwritten if it already exists, then deleted.
        subroutine setup_output()
            character(len=41) :: git_rev

            ! Create NetCDF output file
            call check(nf90_create('output.nc', nf90_clobber, ncid))

            ! Get git revision
            call execute_command_line('git rev-parse HEAD > out.txt')
            open(unit=10, file='out.txt', action='read')
            read(10, *) git_rev
            close(10)
            call execute_command_line('rm -f out.txt')

            ! Write metadata, including model parameters
            call check(nf90_put_att(ncid, nf90_global, "name", "Sam Hatfield"))
            call check(nf90_put_att(ncid, nf90_global, "institution",&
                & "AOPP, University of Oxford"))
            call check(nf90_put_att(ncid, nf90_global, "observation variance", y_var))
            call check(nf90_put_att(ncid, nf90_global, "inflation", rho))
            call check(nf90_put_att(ncid, nf90_global, "localisation", loc))
            call check(nf90_put_att(ncid, nf90_global,&
                & "assimilation frequency",assim_freq))
            call check(nf90_put_att(ncid, nf90_global, "timestep", dt))
            call check(nf90_put_att(ncid, nf90_global, "git-rev", git_rev))

            ! Define time
            call check(nf90_def_dim(ncid, "time", nf90_unlimited, timedim))
            call check(nf90_def_var(ncid, "time", nf90_real4, timedim, timevar))

            ! Define stats variables
            call check(nf90_def_var(ncid, "rmse_mean", nf90_real4, timedim, rmsemeanx))

            ! Write full output?
            if (.not. reduced) then
                call check(nf90_def_dim(ncid, "x", n_x, xdim))
                call check(nf90_def_var(ncid, "x", nf90_int, xdim, xvar))
                call check(nf90_def_dim(ncid, "ens", n_ens, ensdim))
                call check(nf90_def_var(ncid, "ens", nf90_int, ensdim, ensvar))
                call check(nf90_def_var(ncid, "truthx", nf90_real4, (/ xdim, timedim /), truthx))
                call check(nf90_def_var(ncid, "ensx", nf90_real4, (/ xdim, ensdim, timedim /), ensx))
            end if

            call check(nf90_enddef(ncid))
            call check(nf90_close(ncid))
        end subroutine

        subroutine open_file()
            call check(nf90_open('output.nc', nf90_write, ncid))
        end subroutine

        subroutine close_file()
            call check(nf90_close(ncid))
        end subroutine

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Outputs the ensemble and truth data. Only the X values.
        !> @param[in] ensemble the ensemble
        !> @param[in] truth the truth state vector
        !> @param[in] i the time step to write to
        subroutine output(ensemble, truth, i)
            PRECISION, intent(in) :: ensemble(state_dim, n_ens)
            real(dp), intent(in) :: truth(truth_dim)
            integer, intent(in) :: i
            integer :: j

            ! Write stats
            call check(nf90_put_var(ncid, rmsemeanx, rmse_mean(ensemble, truth), (/ i /)))
            call check(nf90_put_var(ncid, timevar, i * dt * write_freq, (/ i /)))
            
            ! Write full output?
            if (.not. reduced) then
                call check(nf90_put_var(ncid, truthx, truth(:n_x), (/ 1, i /)))
                call check(nf90_put_var(ncid, ensx, real(ensemble(:n_x, :)), (/ 1, 1, i /)))
                call check(nf90_put_var(ncid, xvar, (/ (j, j = 1, n_x) /), (/ 1 /)))
                call check(nf90_put_var(ncid, ensvar, (/ (j, j = 1, n_ens) /), (/ 1 /)))
            end if
        end subroutine

        !> @author
        !> Sam Hatfield, AOPP, University of Oxford
        !> @brief
        !> Handles any errors from the NetCDF API.
        !> @param ierr the error code
        subroutine check(ierr)
            integer, intent(in) :: ierr

            if (ierr /= nf90_noerr) then
                print *, trim(adjustl(nf90_strerror(ierr)))
            end if
        end subroutine
end module
