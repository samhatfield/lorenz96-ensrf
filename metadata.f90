module metadata
    use params
    use lorenz96

    implicit none

    contains
        subroutine write_params()
            integer :: file_1 = 20

            open(unit=file_1, file="results.yml", action="write", status="replace")
            1 format(a4, f5.3)
            2 format(a4, i3)
            3 format(a6, i4)
            4 format(a11, i3)
            5 format(a8, f5.2)
            6 format(a4, f6.2)
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
end module metadata
