module local
    use params
    use rp_emulator

    implicit none

    contains
        function localize(gain, y) result(local_gain)
            PRECISION, intent(in) :: gain(obs_dim)
            integer, intent(in) :: y
            PRECISION :: local_gain(obs_dim)
            PRECISION, save :: mask(obs_dim)
            PRECISION :: a
            logical, save :: first_call = .True.
            integer :: i, j

            a = loc

            if (first_call) then
                do i = 0, obs_dim-1
                    ! Account for cyclicity
                    if (i >= obs_dim/2) then
                        j = obs_dim - i
                    else
                        j = i
                    end if
    
                    if (0 <= j .and. j <= a) then
                        mask(i+1) = -0.25_dp*(j/a)**5 + 0.5_dp*(j/a)**4 + 0.625_dp*(j/a)**3 - &
                            & (5.0_dp/3.0_dp)*(j/a)**2 + 1.0_dp
                    else if (a < j .and. j <= 2*a) then
                        mask(i+1) = (1.0_dp/12.0_dp)*(j/a)**5 - 0.5_dp*(j/a)**4 + 0.625_dp*(j/a)**3 + &
                            & (5.0_dp/3.0_dp)*(j/a)**2 - 5.0_dp*(j/a) + 4.0_dp - (2.0_dp/3.0_dp)*(a/j)
                    else
                        mask(i+1) = 0.0_dp
                    end if
                end do

                first_call = .False.
            end if

            local_gain = gain * cshift(mask, 1-y)
        end function
end module
