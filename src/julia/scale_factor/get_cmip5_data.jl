function get_cmip5_data(input_basepath_data)

    # initialize ::
    cmip5_nc, cmip5_data, cmip5_lat, cmip5_lon, cmip5_time_raw, cmip5_time = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()

    # get data ::
    for iage in ["lgm", "piControl"]
        if iage == "lgm"
            idate = "180001-190012"
        elseif iage == "piControl"
            idate = "025001-130012"
        end
        cmip5_nc[iage] = Dict([("u10", NCDataset(joinpath(input_basepath_data, "cmip5/u10/mon_clim/", iage, "/ua_Aclim_CCSM4_" * iage * "_r1i1p1_" * idate * "-clim.nc"))), 
                               ("v10", NCDataset(joinpath(input_basepath_data, "cmip5/u10/mon_clim/", iage, "/va_Aclim_CCSM4_" * iage * "_r1i1p1_" * idate * "-clim.nc")))])
        cmip5_data[iage] = Dict([("u10", core2_nc[iage]["u10"]["U_10"][:]),
                                 ("v10", core2_nc[iage]["v10"]["V_10"][:])])
        cmip5_lat[iage] = Dict([("u10", core2_nc[iage]["u10"]["LAT"][:]),
                                ("v10", core2_nc[iage]["v10"]["LAT"][:])])
        cmip5_lon[iage] = Dict([("u10", core2_nc[iage]["u10"]["LON"][:]),
                                ("v10", core2_nc[iage]["v10"]["LON"][:])])
        cmip5_time_raw[iage] = Dict([("u10", core2_nc[iage]["u10"]["TIME"][:]),
                                     ("v10", core2_nc[iage]["v10"]["TIME"][:])])
        core2_time[iage] = Dict([("u10", ((Dates.dayofyear.(cmip5_time_raw[iage]["u10"]) .- 1) * hr2d) .+ Dates.hour.(cmip5_time_raw[iage]["u10"])),
                                 ("v10", ((Dates.dayofyear.(cmip5_time_raw[iage]["v10"]) .- 1) * hr2d) .+ Dates.hour.(cmip5_time_raw[iage]["v10"]))])
    end

    # return ::
    return cmip5_nc, cmip5_data, cmip5_lat, cmip5_lon, cmip5_time_raw, cmip5_time

end
