module io
    use rp_emulator
    use params
    use netcdf

    implicit none

    private
    public setup_output, output

    integer :: ncid, timedim, xdim, truthx

    contains
        subroutine setup_output()
            call check(nf90_create('results.nc', nf90_clobber, ncid))
            call check(nf90_put_att(ncid, nf90_global, "name", "Sam Hatfield"))
            call check(nf90_put_att(ncid, nf90_global, "institution",&
                & "AOPP, University of Oxford"))
            call check(nf90_def_dim(ncid, "time", nf90_unlimited, timedim))
            call check(nf90_def_dim(ncid, "x", n_x, xdim))
            call check(nf90_def_var(ncid, "truthx", nf90_real8, (/ xdim, timedim /), truthx))
            call check(nf90_enddef(ncid))
            call check(nf90_close(ncid))
        end subroutine

        subroutine output(ensemble, truth)
            PRECISION :: ensemble(state_dim, n_ens)
            real(dp) :: truth(truth_dim)

        end subroutine

        subroutine check(ierr)
            integer, intent(in) :: ierr

            if (ierr /= nf90_noerr) then
                print *, trim(adjustl(nf90_strerror(ierr)))
            end if
        end subroutine
end module
