module metadata
    use params

    implicit none

    contains
        subroutine write_params()
            integer :: file_1 = 20

            open(unit=file_1, file="results.yml", action="write", status="replace")
            1 format(a4, f5.3)
            2 format(a4, i3)
            3 format(a6, i4)
            4 format(a11, i2)
            5 format(a8, f5.2)
            6 format(a4, f6.2)

            write (file_1, 1) 'dt: ', dt 
            write (file_1, 6) 'fin: ', fin
            write (file_1, 2) 'n_x: ', n_x
            write (file_1, 2) 'n_y: ', n_y
            write (file_1, 3) 'n_ens: ', n_ens
            write (file_1, 4) 'assim_freq: ', assim_freq
            write (file_1, 5) 'var_obs: ', var_obs
            write (file_1, '(a10)') '# Columns:'
            write (file_1, '(a42)') '# Max norm, Avg norm, Min norm, Truth, Obs'
            write (file_1, '(a5)') '--- |'

            close(file_1)
        end subroutine write_params
end module metadata