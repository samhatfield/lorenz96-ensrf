module observation
    use rp_emulator
    use params

    implicit none

    !===========================================================================
    ! Overloaded functions
    !=========================================================================== 

    public :: observe
    interface observe
        module procedure observe
        module procedure observe_rpe
        module procedure observe_row
        module procedure observe_row_rpe
    end interface observe

    interface observe_1d_row
        module procedure observe_1d_row
        module procedure observe_1d_row_rpe
    end interface observe_1d_row

    !===========================================================================
    ! Function definitions
    !=========================================================================== 

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

        pure function observe_row(state, row)
            real(dp), dimension(:, :), intent(in) :: state
            integer, intent(in) :: row
            real(dp), dimension(size(state, 2)) :: observe_row

            observe_row = state(n_x+row,:)
        end function observe_row

        pure function observe_row_rpe(state, row)
            type(rpe_var), dimension(:, :), intent(in) :: state
            integer, intent(in) :: row
            type(rpe_var), dimension(size(state, 2)) :: observe_row_rpe

            observe_row_rpe = observe_row(state%val, row)
        end function observe_row_rpe

        pure function observe_1d_row(state, row)
            real(dp), dimension(state_dim), intent(in) :: state
            integer, intent(in) :: row
            real(dp) :: observe_1d_row

            observe_1d_row = state(n_x+row)
        end function observe_1d_row

        pure function observe_1d_row_rpe(state, row)
            type(rpe_var), dimension(state_dim), intent(in) :: state
            integer, intent(in) :: row
            type(rpe_var) :: observe_1d_row_rpe

            observe_1d_row_rpe = observe_1d_row(state%val, row)
        end function observe_1d_row_rpe
end module observation
