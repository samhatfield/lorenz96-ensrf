include("params.jl")

const F = 20.0
const h = 1.0
const c = 10.0
const b = 10.0
const e = 10.0
const d = 10.0
const gZ = 1.0

"Step forward once using 4th order Runge-Kutta"
function step(r, rhs)
    k₁ = rhs(r)
    k₂ = rhs(r + 0.5Δt*k₁)
    k₃ = rhs(r + 0.5Δt*k₂)
    k₄ = rhs(r + Δt*k₃)

    r + (Δt/6.0)*(k₁ + 2k₂ + 2k₃ + k₄)
end

"Definition of RHS for Lorenz three-level '96 system'"
function three_level_rhs(r)
    # Split up state vector into components
    x = r[1:nx]
    y = r[nx+1:nx+nx*ny]
    z = r[nx+nx*ny+1:end]

    # Find derivative of each component separately
    vcat(dxdt(x, y), dydt(x, y, z), dzdt(y, z))
end

"Definition of RHS for Lorenz two-level '96 system"
function two_level_rhs(r)
    # Split up state vector into components
    x = r[1:nx]
    y = r[nx+1:end]

    # Find derivative of each component separately
    vcat(dxdt(x, y), dydt(x, y))
end

"Definition of RHS for X variables in three-level Lorenz '96 system"
function dxdt(x, y)
    sum_y = dropdims(sum(reshape(y, (ny,nx)), dims=1), dims=1)

    circshift(x,-1).*(circshift(x,1) - circshift(x,-2)) - x - (h*c/b)*sum_y .+ F
end

"Definition of RHS for Y variables in three-level Lorenz '96 system"
function dydt(x, y, z)
    # Repeat elements of x ny times
    x_rpt = [x[1+floor(Int, (i-1)/ny)] for i=1:nx*ny]

    # Sum all z's for each y, making an nx*ny length vector of z sums
    sum_z = dropdims(sum(reshape(z, (nz,nx*ny)), dims=1), dims=1)

    c*b*circshift(y,1).*(circshift(y,-1) - circshift(y,2)) - c*y + (h*c/b)*x_rpt - (h*e/d)*sum_z
end

"Definition of RHS for Y variables in two-level Lorenz '96 system"
function dydt(x, y)
    # Repeat elements of x ny times
    x_rpt = [x[1+floor(Int, (i-1)/ny)] for i=1:nx*ny]

    c*b*circshift(y,1).*(circshift(y,-1) - circshift(y,2)) - c*y + (h*c/b)*x_rpt
end

"Definition of RHS for Z variables in three-level Lorenz '96 system'"
function dzdt(y, z)
    # Repeat elements of y nz times
    y_rpt = [y[1+floor(Int, (i-1)/nz)] for i=1:nx*ny*nz]

    e*d*circshift(z,-1).*(circshift(z,1) - circshift(z,-2)) - gZ*e*z + (h*e/d)*y_rpt
end
