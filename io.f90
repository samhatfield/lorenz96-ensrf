!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains functions necessary for writing to an output NetCDF file.
module io
    use rp_emulator
    use params
    use netcdf
    use utils, only: real, rmse_ens_mean, mean_ens_rmse

    implicit none

    private
    public setup_output, output, open_file, close_file

    integer :: ncid
    integer :: timedim, xdim, ensdim, ydim
    integer :: timevar, xvar, truthx, ensvar, ensx, yvar, truthy, ensy, rmseensmeanx, meanensrmsex

    contains
        !> @brief
        !> Sets up the output NetCDF file, including writing simulation
        !> parameters.
        subroutine setup_output()
            ! Create NetCDF output file
            call check(nf90_create(trim(outfile), nf90_clobber, ncid))

            ! Write metadata, including model parameters
            call check(nf90_put_att(ncid, nf90_global, "name", "Sam Hatfield"))
            call check(nf90_put_att(ncid, nf90_global, "institution",&
                & "AOPP, University of Oxford"))
            call check(nf90_put_att(ncid, nf90_global, "timestep", dt))
            call check(nf90_put_att(ncid, nf90_global, "x variable count", n_x))
            call check(nf90_put_att(ncid, nf90_global, "y variable count (per x)", n_y))
            call check(nf90_put_att(ncid, nf90_global, "z variable count (per y)", n_z))
            call check(nf90_put_att(ncid, nf90_global, "ensemble size", n_ens))
            call check(nf90_put_att(ncid, nf90_global,&
                & "assimilation frequency",assim_freq))
            call check(nf90_put_att(ncid, nf90_global, "observation variance", y_var))
            call check(nf90_put_att(ncid, nf90_global, "observation dimensions", obs_dim))
            call check(nf90_put_att(ncid, nf90_global, "inflation", rho))
            call check(nf90_put_att(ncid, nf90_global, "localisation", loc))
            call check(nf90_put_att(ncid, nf90_global, "main git rev", GIT_REV_MAIN))
            call check(nf90_put_att(ncid, nf90_global, "rpe git rev", GIT_REV_RPE))
            call check(nf90_put_att(ncid, nf90_global, "precision", PREC_STR))

            ! Define time
            call check(nf90_def_dim(ncid, "time", nf90_unlimited, timedim))
            call check(nf90_def_var(ncid, "time", nf90_real4, timedim, timevar))

            ! Define stats variables
            call check(nf90_def_var(ncid, "rmse_ens_mean", nf90_real4, timedim, rmseensmeanx))
            call check(nf90_def_var(ncid, "mean_ens_rmse", nf90_real4, timedim, meanensrmsex))

            ! Output X or Y?
            if (output_x .or. output_y) then
                call check(nf90_def_dim(ncid, "ens", n_ens, ensdim))
                call check(nf90_def_var(ncid, "ens", nf90_int, ensdim, ensvar))
            end if

            ! Output X?
            if (output_x) then
                call check(nf90_def_dim(ncid, "x", n_x, xdim))
                call check(nf90_def_var(ncid, "x", nf90_int, xdim, xvar))
                call check(nf90_def_var(ncid, "truthx", nf90_real4, (/ xdim, timedim /), truthx))
                call check(nf90_def_var(ncid, "ensx", nf90_real4, (/ xdim, ensdim, timedim /), ensx))
            end if

            ! Write X and Y
            if (output_y) then
                call check(nf90_def_dim(ncid, "y", n_x*n_y, ydim))
                call check(nf90_def_var(ncid, "y", nf90_int, ydim, yvar))
                call check(nf90_def_var(ncid, "truthy", nf90_real4, (/ ydim, timedim /), truthy))
                call check(nf90_def_var(ncid, "ensy", nf90_real4, (/ ydim, ensdim, timedim /), ensy))
            end if

            call check(nf90_enddef(ncid))
            call check(nf90_close(ncid))
        end subroutine

        subroutine open_file()
            call check(nf90_open(outfile, nf90_write, ncid))
        end subroutine

        subroutine close_file()
            real(dp) :: elapsed

            ! Get elapsed time
            call cpu_time(elapsed)
            call check(nf90_redef(ncid))
            call check(nf90_put_att(ncid, nf90_global, "cpu_time_seconds", elapsed))

            call check(nf90_close(ncid))
        end subroutine

        !> @brief
        !> Outputs the ensemble and truth data. Only the X values.
        !> @param[in] ensemble the ensemble
        !> @param[in] truth the truth state vector
        !> @param[in] i the time step to write to
        subroutine output(ensemble, truth, i)
            PRECISION, intent(in) :: ensemble(state_dim, n_ens)
            real(dp), intent(in) :: truth(state_dim)
            integer, intent(in) :: i
            integer :: j

            ! Write stats
            call check(nf90_put_var(ncid, timevar, i * dt * write_freq, (/ i /)))
            call check(nf90_put_var(ncid, rmseensmeanx, rmse_ens_mean(ensemble, truth), (/ i /)))
            call check(nf90_put_var(ncid, meanensrmsex, mean_ens_rmse(ensemble, truth), (/ i /)))
            
            ! Output X or Y?
            if (output_x .or. output_y) then
                call check(nf90_put_var(ncid, ensvar, (/ (j, j = 1, n_ens) /), (/ 1 /)))
            end if

            ! Output X?
            if (output_x) then
                call check(nf90_put_var(ncid, truthx, truth(:n_x), (/ 1, i /)))
                call check(nf90_put_var(ncid, ensx, real(ensemble(:n_x, :)), (/ 1, 1, i /)))
                call check(nf90_put_var(ncid, xvar, (/ (j, j = 1, n_x) /), (/ 1 /)))
            end if

            ! Output Y?
            if (output_y) then
                call check(nf90_put_var(ncid, truthy, truth(n_x+1:n_x+n_x*n_y), (/ 1, i /)))
                call check(nf90_put_var(ncid, ensy, real(ensemble(n_x+1:n_x+n_x*n_y, :)), (/ 1, 1, i /)))
                call check(nf90_put_var(ncid, yvar, (/ (j, j = 1, n_x*n_y) /), (/ 1 /)))
            end if
        end subroutine

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
