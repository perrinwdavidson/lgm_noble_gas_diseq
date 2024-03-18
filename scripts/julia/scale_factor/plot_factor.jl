function plot_factor(core2_lat, cmip5_zonal_scale_factor, output_basepath_plots)

    i, j = 1, 1
    f = Figure(resolution = (1800, 1200))
    for imonth in 1 : 1 : size(cmip5_data["lgm"][ivar], 3)
        scatterlines(f[i, j], core2_lat[ivar], cmip5_zonal_scale_factor["u10"][:, imonth], axis = (; title = "Month " * string(imonth), xlabel = "Latitude (deg. North)", ylabel = "Wind Factor"))
        j += 1
        if j > 3
            j = 1
            i += 1
        end
    end
    save(output_basepath_plots * "/u10/scale_factor/scale_factor.png", f, px_per_unit = 10)

end
