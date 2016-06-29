module observation
    use rp_emulator
    use params

    implicit none

    !===========================================================================
    ! Overloaded functions
    !=========================================================================== 

    public :: observe
    interface observe
        {% for type in types %}
        module procedure observe_{{ type.name }}
        module procedure observe_row_{{ type.name }}
        module procedure observe_1d_row_{{ type.name }}
        {% endfor %}
    end interface

    !===========================================================================
    ! Function definitions
    !=========================================================================== 

    contains
        {% for type in types %}
        pure function observe_{{ type.name }}(state) result(observe)
            {{ type.code }}, intent(in) :: state(:,:)
            {{ type.code }} :: observe(obs_dim, size(state,2))
            
            observe = state(n_x+1:n_x+n_x*n_y,:)
        end function

        pure function observe_row_{{ type.name }}(state, row) result(observe_row)
            {{ type.code }}, intent(in) :: state(:,:)
            integer, intent(in) :: row
            {{ type.code }} :: observe_row(size(state,2))

            observe_row = state(n_x+row,:)
        end function

        pure function observe_1d_row_{{ type.name }}(state, row) result(observe_1d_row)
            {{ type.code }}, intent(in) :: state(state_dim)
            integer, intent(in) :: row
            {{ type.code }} :: observe_1d_row

            observe_1d_row = state(n_x+row)
        end function
        {% endfor %}
end module
