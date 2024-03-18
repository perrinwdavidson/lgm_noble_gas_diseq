# =========================================================
# interp_pmip3
# ---------------------------------------------------------
# purpose :: to interpolate pmip3 data to uvic grid
# author :: perrin w. davidson
# contact :: perrinwdavidson@gmail.com
# date :: 05.04.23
# =========================================================
# configure -----------------------------------------------
# load developmental geokrige ::
include("/Users/perrindavidson/.julia/dev/GeoKrige/src/GeoKrige.jl");  # remove after successful development 

# load packages ::
using .GeoKrige  # for kriging
using NCDatasets  # for netcdf files
using DelimitedFiles  # for loading delimited files
using Statistics  # for statistical functionality
using PyPlot  # for plotting
using PyCall  # for calling python
using Random  # for random functionality
using GroupSlices  # for array functionality
using ScatteredInterpolation  # for interpolating functionality

# set working directory ::
cd("/Users/perrindavidson/research/whoi/current/noble_gas"); 

# specify plotting ::
plt.rc("text", usetex=true);
plt.rc("font", family="serif", size=16);
axes_grid1 = pyimport("mpl_toolkits.axes_grid1");
colors = pyimport("matplotlib.colors"); 

# set basepath ::
basepath = "/Users/perrindavidson/Library/CloudStorage/GoogleDrive-perrinwd@mit.edu/My Drive/research/whoi/current/noble_gas/cmip5/";
basepath_query = joinpath(basepath, "data/exp_raw/uvic/grid/");

# interpolate data ----------------------------------------
# set age types ::
age_types = ["pic", "lgm"]; 

# set variable types ::
variable_types = ["sic", "u10"];

# loop over all data ::
for iage in age_types
	for ivar in variable_types

		# set basepaths :
		basepath_sample = joinpath(basepath, "data", "exp_raw", "cmip5", ivar, "mon_clim", iage);
		basepath_plots = joinpath(basepath, "plots", ivar, iage);
		basepath_interp = joinpath(basepath, "data", "sims", ivar, iage);
	
		# get filenames ::
		filenames = readdir(basepath_sample);
	
		# determine interpolation function ::
		if ivar == "sic"
			for fname in filenames
				if !isnothing(findfirst("sic", fname)) && isnothing(findfirst("get", fname))

					# get model name ::
					model_name = lowercase(split(split(fname, "_")[3], "-")[1]);					

					# interpolate ::
					eval(Meta.parse(replace("include('" * model_name * "/interp_sic_" * model_name * ".jl')", "'" => '"')));
					eval(Meta.parse("interp_sic_" * model_name * "(fname, basepath_sample, basepath_plots, basepath_interp)"));

				else
					continue
				end
			end
		elseif ivar == "u10"
			for fname in filenames
				if (!isnothing(findfirst("ua", fname)) || !isnothing(findfirst("va", fname))) && isnothing(findfirst("get", fname))
					interp_u10(fname, basepath_sample, basepath_plots, basepath_interp);
				else
					continue
				end
			end
		end

	end
end

# =========================================================
