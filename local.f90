module local
    use params
    use rp_emulator

    implicit none

    contains
        function localize_y(gain, y) result(local_gain)
            PRECISION, intent(in) :: gain(obs_dim)
            integer, intent(in) :: y
            PRECISION :: local_gain(obs_dim)
            PRECISION, save :: mask(obs_dim)
            PRECISION :: c
            logical, save :: first_call = .True.
            integer :: i, j

            c = loc_y

            if (first_call) then
                do i = 0, obs_dim-1
                    ! Account for cyclicity
                    if (i >= obs_dim/2) then
                        j = obs_dim - i
                    else
                        j = i
                    end if
    
                    mask(i+1) = gaspari_cohn(j, c)
                end do

                first_call = .False.
            end if

            local_gain = gain * cshift(mask, 1-y)
        end function

        function localize_x(gain, y) result(local_gain)
            PRECISION, intent(in) :: gain(n_x)
            integer, intent(in) :: y
            PRECISION :: local_gain(n_x)
            PRECISION, save :: mask(n_x)
            PRECISION :: c
            logical, save :: first_call = .True.
            integer :: i, j
            integer :: dist_y

            dist_y = int((y-1)/n_y)+1

            c = loc_x

            if (first_call) then
                do i = 0, n_x-1
                    ! Account for cyclicity
                    if (i >= n_x/2) then
                        j = n_x - i
                    else
                        j = i
                    end if

                    mask(i+1) = gaspari_cohn(j, c)
                end do

                first_call = .False.
            end if

            local_gain = gain * cshift(mask, 1-dist_y)
        end function

        function gaspari_cohn(z_, c) result(C_0)
            integer, intent(in) :: z_
            PRECISION, intent(in) :: c
            PRECISION :: z
            PRECISION :: C_0

            z = z_

            if (0 <= z .and. z <= c) then
                C_0 = -0.25_dp*(z/c)**5 + 0.5_dp*(z/c)**4 + 0.625_dp*(z/c)**3 - &
                    & (5.0_dp/3.0_dp)*(z/c)**2 + 1.0_dp
            else if (c < z .and. z <= 2*c) then
                C_0 = (1.0_dp/12.0_dp)*(z/c)**5 - 0.5_dp*(z/c)**4 + 0.625_dp*(z/c)**3 + &
                    & (5.0_dp/3.0_dp)*(z/c)**2 - 5.0_dp*(z/c) + 4.0_dp - (2.0_dp/3.0_dp)*(c/z)
            else
                C_0 = 0.0_dp
            end if
        end function
end module
