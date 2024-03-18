# =========================================================
# calc_scale_factor
# ---------------------------------------------------------
# purpose :: to calculate the scale factor for core-2 winds
# author :: perrin w. davidson
# contact :: perrinwdavidson@gmail.com
# date :: 12.06.23
# =========================================================
# configure -----------------------------------------------
# load packages ::
using Pkg
using Revise
using NCDatasets
using CFTime
using Dates
using Statistics
using Interpolations
using CairoMakie

# include functions ::
include("get_cmip5_data.jl")
include("get_core2_data.jl")
include("plot_factor.jl")

# environment ::
cd("/Users/perrindavidson/research/whoi/current/noble_gas")

# set basepaths ::
basepath = "/Users/perrindavidson/Library/CloudStorage/GoogleDrive-perrinwd@mit.edu/My Drive/research/whoi/current/noble_gas/cmip5"
input_basepath_data = joinpath(basepath, "data/exp_raw")
output_basepath_data = joinpath(basepath, "data/sims")
output_basepath_plots = joinpath(basepath, "plots")

# load data -----------------------------------------------
# core-2 data ::
core2_nc, core2_data, core2_lat, core2_lon, core2_time_raw, core2_time = get_core2_data(input_basepath_data)

# cmip5 data ::
cmip5_nc, cmip5_data, cmip5_lat, cmip5_lon, cmip5_time_raw, cmip5_time = get_cmip5_data(input_basepath_data)

# make factor ---------------------------------------------
# calculate cumulative sum ::
cmip5_zonal_scale_factor_data = Dict()
for ivar in ["u10", "v10"]
    cmip5_zonal_scale_factor_data[ivar] = Dict()
    for iage in ["lgm", "piControl"]
        cmip5_zonal_scale_factor_data[ivar][iage] = Dict([("sum", zeros(Float64, length(cmip5_lat[iage][ivar]))), ("n", [])])
        for imonth = 1 : 1 : size(cmip5_data[iage][ivar], 3)
            cmip5_zonal_scale_factor_data[ivar][iage]["sum"][imonth] = sum.(skipmissing.(eachrow(cmip5_data[iage][ivar][:, :, imonth]')))
            cmip5_zonal_scale_factor_data[ivar][iage]["n"] = length(skipmissing.(eachrow(cmip5_data[iage][ivar][:, :, imonth]')))
        end        
    end
end

# plot put ::
# plot_factor(core2_lat, cmip5_zonal_scale_factor, output_basepath_plots)

 
#cmip5_zonal_scale_factor = Dict()
#for ivar in ["u10", "v10"]
#    cmip5_zonal_scale_factor[ivar] = zeros(Float64, length(core2_lat[ivar]), size(cmip5_data["lgm"][ivar], 3))
#    for imonth = 1 : 1 : size(cmip5_data["lgm"][ivar], 3)
#        dat = mean.(skipmissing.(eachrow(cmip5_data["lgm"][ivar][:, :, imonth]'))) ./ mean.(skipmissing.(eachrow(cmip5_data["piControl"][ivar][:, :, imonth]')))
#        lin_interp = linear_interpolation(cmip5_lat["lgm"][ivar], dat)
#        cmip5_zonal_scale_factor[ivar][:, imonth] = lin_interp.(core2_lat[ivar])
#    end        
#end


# =========================================================
