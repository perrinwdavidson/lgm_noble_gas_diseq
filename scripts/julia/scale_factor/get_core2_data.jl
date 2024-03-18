function get_core2_data(input_basepath_data)

    # set time constant ::
    hr2d = 24  # [hours in a day]
    d2yr = 365  # [days in a year]
    hr2yr = hr2d * d2yr  # [hours in a year]

    # initialize ::
    core2_nc, core2_data, core2_lat, core2_lon, core2_time_raw, core2_time = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
    
    # get data ::
    core2_nc = Dict([("u10", NCDataset(joinpath(input_basepath_data, "core2/u_10.15JUNE2009.nc"))), 
                     ("v10", NCDataset(joinpath(input_basepath_data, "core2/v_10.15JUNE2009.nc")))])
    core2_data = Dict([("u10", core2_nc["u10"]["U_10"][:]),
                       ("v10", core2_nc["v10"]["V_10"][:])])
    core2_lat = Dict([("u10", core2_nc["u10"]["LAT"][:]),
                      ("v10", core2_nc["v10"]["LAT"][:])])
    core2_lon = Dict([("u10", core2_nc["u10"]["LON"][:]),
                      ("v10", core2_nc["v10"]["LON"][:])])
    core2_time_raw = Dict([("u10", core2_nc["u10"]["TIME"][:]),
                           ("v10", core2_nc["v10"]["TIME"][:])])
    core2_time = Dict([("u10", ((Dates.dayofyear.(core2_time_raw["u10"]) .- 1) * hr2d) .+ Dates.hour.(core2_time_raw["u10"])),
                       ("v10", ((Dates.dayofyear.(core2_time_raw["v10"]) .- 1) * hr2d) .+ Dates.hour.(core2_time_raw["v10"]))])

    # return ::
    return core2_nc, core2_data, core2_lat, core2_lon, core2_time_raw, core2_time

end
