!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> Contains observation operator in several forms e.g. for generating a whole
!> observation vector, or just a single scalar observation.
module observation
    use rp_emulator
    use params

    implicit none

    private

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
        !> @brief
        !> Observe a whole ensemble of state vectors. This treats each column
        !> (or row, I forget) as a state vector, observes each and returns
        !> a matrix of observed state vectors.
        !> @param state the state vector
        !> @return observe the observed matrix
        pure function observe_{{ type.name }}(state) result(observe)
            {{ type.code }}, intent(in) :: state(:,:)
            {{ type.code }} :: observe(obs_dim, size(state,2))
            
            observe = state(n_x+1:n_x+n_x*n_y:y_skip,:)
        end function

        !> @brief
        !> Observe a whole ensemble of state vectors, but just a single
        !> Y variable, instead of all observed variables. This treats each
        !> column (or row, I forget) as a state vector, observes the Y given
        !> variable in each and returns a vector of observed variables drawn
        !> from each state vector. This is used for sequential data
        !> assimilation, where you process only one observed variable at a time.
        !> @param state the state vector
        !> @param y the Y variable to be observed
        !> @return observe the observed matrix
        pure function observe_row_{{ type.name }}(state, y) result(observe_row)
            {{ type.code }}, intent(in) :: state(:,:)
            integer, intent(in) :: y
            {{ type.code }} :: observe_row(size(state,2))

            observe_row = state(n_x+(y_skip*y-(y_skip-1)),:)
        end function

        !> @brief
        !> Observe a single state vector, but just a single Y variable, instead
        !> of all observed variables. This is used for sequential data
        !> assimilation, where you process only one observed variable at a time.
        !> @param state the state vector
        !> @param y the Y variable to be observed
        !> @return observe the observed matrix
        pure function observe_1d_row_{{ type.name }}(state, y) result(observe_1d_row)
            {{ type.code }}, intent(in) :: state(state_dim)
            integer, intent(in) :: y 
            {{ type.code }} :: observe_1d_row

            observe_1d_row = state(n_x+(y_skip*y-(y_skip-1)))
        end function
        {% endfor %}
end module
