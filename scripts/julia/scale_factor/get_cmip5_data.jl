function get_cmip5_data(input_basepath_data)

    # set time constant ::
    hr2d = 24  # [hours in a day]
    d2yr = 365  # [days in a year]
    hr2yr = hr2d * d2yr  # [hours in a year]

    # initialize ::
    cmip5_nc, cmip5_data, cmip5_lat, cmip5_lon, cmip5_time_raw, cmip5_time = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()

    # get data ::
    for iage in ["lgm", "piControl"]
        if iage == "lgm"
            idate = "180001-190012"
        elseif iage == "piControl"
            idate = "025001-130012"
        end
        fname_u10 = "cmip5/u10/mon_clim/" * iage * "/ua_Aclim_CCSM4_" * iage * "_r1i1p1_" * idate * "-clim.nc"
        fname_v10 = "cmip5/u10/mon_clim/" * iage * "/va_Aclim_CCSM4_" * iage * "_r1i1p1_" * idate * "-clim.nc"
        cmip5_nc[iage] = Dict([("u10", NCDataset(joinpath(input_basepath_data, fname_u10))), 
                               ("v10", NCDataset(joinpath(input_basepath_data, fname_v10)))])
        cmip5_data[iage] = Dict([("u10", cmip5_nc[iage]["u10"]["ua"][:, :, 1, :]),  # [100000.0 Pa]
                                 ("v10", cmip5_nc[iage]["v10"]["va"][:, :, 1, :])])  # [100000.0 Pa]
        cmip5_lat[iage] = Dict([("u10", cmip5_nc[iage]["u10"]["lat"][:]),
                                ("v10", cmip5_nc[iage]["v10"]["lat"][:])])
        cmip5_lon[iage] = Dict([("u10", cmip5_nc[iage]["u10"]["lon"][:]),
                                ("v10", cmip5_nc[iage]["v10"]["lon"][:])])
        cmip5_time_raw[iage] = Dict([("u10", cmip5_nc[iage]["u10"]["time"][:]),
                                     ("v10", cmip5_nc[iage]["v10"]["time"][:])])
        cmip5_time[iage] = Dict([("u10", ((Dates.dayofyear.(cmip5_time_raw[iage]["u10"]) .- 1) * hr2d) .+ Dates.hour.(cmip5_time_raw[iage]["u10"])),
                                 ("v10", ((Dates.dayofyear.(cmip5_time_raw[iage]["v10"]) .- 1) * hr2d) .+ Dates.hour.(cmip5_time_raw[iage]["v10"]))])
    end

    # return ::
    return cmip5_nc, cmip5_data, cmip5_lat, cmip5_lon, cmip5_time_raw, cmip5_time

end
