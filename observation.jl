#= Observe a whole ensemble of state vectors. This treats each column
   (or row, I forget) as a state vector, observes each and returns
   a matrix of observed state vectors. =#
H(state) = state[nx+1:nx+nx*ny,:]

#= Observe a whole ensemble of state vectors, but just a single
   Y variable, instead of all observed variables. This treats each
   column (or row, I forget) as a state vector, observes the Y given
   variable in each and returns a vector of observed variables drawn
   from each state vector. This is used for sequential data
   assimilation, where you process only one observed variable at a time. =#
function H(state, row)
    # If observing an ensemble
    if length(size(state)) == 2
        return state[nx+row,:]
    # Else if observing just one state vector
    else
        return state[nx+row]
    end
end
