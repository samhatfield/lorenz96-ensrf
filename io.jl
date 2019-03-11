include("params.jl")

using NetCDF

"Sets up the output NetCDF file, including writing simulation parameters."
function set_up_output()
    # Define dimensions
    timedim = NcDim("time", 0, unlimited=true)
    xdim    = NcDim("x", nx)
    ydim    = NcDim("y", nx*ny)
    ensdim  = NcDim("ens", nens)

    # Define dimension variables
    timevar = NcVar("time", timedim, t=Int32)
    xvar    = NcVar("x",    xdim,    t=Int32)
    yvar    = NcVar("y",    ydim,    t=Int32)
    ensvar  = NcVar("ens",  ensdim,  t=Int32)

    # Define multidimensional variables
    truthx = NcVar("truthx",   [xdim, timedim], t=Float32)
    truthy = NcVar("truthy",   [ydim, timedim], t=Float32)
    ensx   = NcVar("ensx",     [xdim, ensdim, timedim])
    ensy   = NcVar("ensy",     [ydim, ensdim, timedim])

    nc = NetCDF.create(outfile, [truthx, truthy, ensx, ensy, timevar], mode=NC_NETCDF4)
end

"Outputs the ensemble and truth data."
function output(ncfile, t_index, ensemble, truth)
    NetCDF.putvar(ncfile, "truthx", truth[1:nx],          start=[1,t_index], count=[-1,1])
    NetCDF.putvar(ncfile, "truthy", truth[nx+1:nx+nx*ny], start=[1,t_index], count=[-1,1])
    NetCDF.putvar(ncfile, "ensx", ensemble[1:nx,:],          start=[1,1,t_index], count=[-1,-1,1])
    NetCDF.putvar(ncfile, "ensy", ensemble[nx+1:nx+nx*ny,:], start=[1,1,t_index], count=[-1,-1,1])
end
