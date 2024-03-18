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

# specify plotting ::
plt.rc("text", usetex=true);
plt.rc("font", family="serif", size=16);
axes_grid1 = pyimport("mpl_toolkits.axes_grid1");
colors = pyimport("matplotlib.colors"); 

# load data -----------------------------------------------
# basepath ::
basepath_sample = "/Users/perrindavidson/Library/CloudStorage/GoogleDrive-perrinwd@mit.edu/My Drive/research/whoi/current/noble_gas/cmip5/data/exp_raw/cmip5/sic/mon_clim/pic/";
basepath_query = "/Users/perrindavidson/Library/CloudStorage/GoogleDrive-perrinwd@mit.edu/My Drive/research/whoi/current/noble_gas/cmip5/data/exp_raw/uvic/grid/";
basepath_plots = "/Users/perrindavidson/Library/CloudStorage/GoogleDrive-perrinwd@mit.edu/My Drive/research/whoi/current/noble_gas/cmip5/plots/sic/pic/";
basepath_krige = "/Users/perrindavidson/Library/CloudStorage/GoogleDrive-perrinwd@mit.edu/My Drive/research/whoi/current/noble_gas/cmip5/data/sims/sic/pic/";

# sample data ::
fname = "sic_OIclim_CCSM4_piControl_r1i1p1_125001-230012-clim.nc";
sample_nc = NCDataset(joinpath(basepath_sample, fname), "r");
sample_x_grid = Matrix{Float64}(sample_nc["lon"][:]);
sample_y_grid = Matrix{Float64}(sample_nc["lat"][:]);
sample_f_grid = sample_nc["sic"][:, :, :];
sample_x, sample_y = sample_x_grid[:], sample_y_grid[:];

# get model name ::
model_name = split(fname, "_")[3]; 

# make longitude [0, 360] ::
sample_x = shift_lon(sample_x, sample_x); 

# query grid ::
query_x_grid = readdlm(joinpath(basepath_query, "grid_lon.csv"), ',', Float64);
query_y_grid = readdlm(joinpath(basepath_query, "grid_lat.csv"), ',', Float64);
query_mask_grid = BitArray(readdlm(joinpath(basepath_query, "grid_mask.csv"), ',', Float64)); 

# get number of months ::
num_months = size(sample_f_grid)[3];

# set ice bound ::
ice_lat_bound = -65;

# make directory if it doesn't exist ::
if !isdir(joinpath(basepath_plots, model_name))
	mkdir(joinpath(basepath_plots, model_name)); 
end
if !isdir(joinpath(basepath_krige, model_name))
	mkdir(joinpath(basepath_krige, model_name)); 
end

# preallocate final array ::
pmip3_data = zeros(size(query_mask_grid)[1], size(query_mask_grid)[2], num_months); 
pmip3_var = zeros(size(query_mask_grid)[1], size(query_mask_grid)[2], num_months); 

# print out ::
println("Done loading data.");

# loop through all months ---------------------------------
for imonth in 1 : 1 : num_months

	# quality control data ----------------------------
	# get monthly data ::
	sample_f = sample_f_grid[:, :, imonth];

	# flatten query grid ::
	query_x = query_x_grid[:]; 
	query_y = query_y_grid[:]; 
	query_mask = query_mask_grid[:]; 
	query = hcat(query_x, query_y); 

	# get nan idx ::
	idx_sample_real = findall(!ismissing, sample_f); 

	# make in [0, 1] domain ::
	if mean(sample_f[idx_sample_real]) > 1
		sample_f ./= 100;
	end

	# artificially impose ice -------------------------
	# find missing ::
	idx_nan = findall(ismissing, sample_f);

	# find above lat degrees ::
	idx_lat_bound = findall(<=(ice_lat_bound), sample_y_grid);  

	# find intersection of both ::
	idx_ice = intersect(idx_lat_bound, idx_nan); 

	# impose ice ::
	sample_f[idx_nan] .= NaN;
	sample_f[idx_ice] .= 1.;

	# update indices ::
	idx_sample_real = findall(!isnan, sample_f[:]); 

	# plot ::
	# fig, ax = plt.subplots(1, 1, figsize=[12, 6]); 
	# cax = ax.scatter(sample_x_grid, sample_y_grid, c=sample_f); 
	# plt.colorbar(cax); 
	# plt.tight_layout();
	# plt.show();

	# make inputs -------------------------------------
	# normalize and detrend data ::
	sample_f_norm_detrend, sample_f_norm_mean = detrend(sample_f[idx_sample_real]); 
	
	# make sample array ::
	sample = hcat(sample_x[idx_sample_real], sample_y[idx_sample_real], sample_f_norm_detrend); 

	# plot normalized data ::
	fig, ax = plt.subplots(1, 2, figsize=[12, 6]);
	ax[1].hist(sample_f, bins=sturges(length(sample_f))); 
	ax[1].set_title("Raw");
	ax[1].set_ylabel("Counts");
	ax[1].set_xlabel("Data");
	ax[2].hist(sample_f_norm_detrend, bins=sturges(length(sample_f_norm_detrend))); 
	ax[2].set_title("Normalized");
	ax[2].set_ylabel("Counts");
	ax[2].set_xlabel("Data");
	plt.tight_layout(); 
	plt.savefig(joinpath(basepath_plots, model_name, model_name * "_" * string(imonth) * "_norm_detrend_data.png"), dpi=300); 
	plt.show(); 

	# print out ::
	println("Done making inputs.");

	# calculate experimental variogram ----------------
	# set inputs ::
	max_lag = 5E6; 
	lags, semivar = Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}();

	# choose subset of data ::
	idx_subset = randsubseq(MersenneTwister(1234), 1 : 1 : length(sample_f_norm_detrend), 0.1)

	# calculate matheron estimator ::
	ilags, isemivar = semivariogram(sample[idx_subset, 1:2], 
					sample[idx_subset, 3], 
					max_lag; 
					estimator="matheron");
	lags = push!(lags, "matheron" => ilags); 
	semivar = push!(semivar, "matheron" => isemivar); 

	# calculate cressie estimator ::
	# ilags, isemivar = semivariogram(sample[idx_subset, 1:2], 
	#				  sample[idx_subset, 3], 
	#				  max_lag; 
	#				  estimator="cressie");
	# lags = push!(lags, "cressie" => ilags); 
	# semivar = push!(semivar, "cressie" => isemivar); 

	# plot ::
	# fig, ax = plt.subplots(1, 1, figsize=[6, 6]);
	# ax.scatter(lags["matheron"], semivar["matheron"], label="Matheron");
	# ax.scatter(lags["cressie"], semivar["cressie"], label="Cressie");
	# ax.set_xlabel("Lag [m]");
	# ax.set_ylabel("Semivariance"); 
	# ax.set_title("Semivariogram");
	# plt.legend();
	# plt.tight_layout(); 
	# plt.show();

	# fit semivariogram -------------------------------
	# make inputs ::
	vario_estimator = "matheron"; 
	p0 = [0., 1E-2, 2E6];

	# choose models ::
	models =["Spherical", "Gaussian", "Exponential", "Power"]; 

	# fit variogram :: 
	best_model, best_model_name, vario_err, pf, pferr, chisq, dof = fit_variogram(lags[vario_estimator], 
										      semivar[vario_estimator], 
										      p0, 
										      models; 
										      plot_output=true,
										      save_fig=joinpath(basepath_plots, model_name, model_name * "_" * string(imonth) * "_fit_semivariogram.png"));		

	# display ::
	println("Done calculating and fitting experimental variogram.");
	println("Best model is: ", best_model_name); 

	# krige -------------------------------------------
	# set parameters ::
	max_lag_krige = 5E5; 

	# krige data ::
	krige_data_norm_detrend = krige(sample,
					query,
					query_mask,
					best_model,
					pf[best_model_name],
					max_lag_krige); 	

	# renormalize and add in mean ::
	krige_data = deepcopy(krige_data_norm_detrend); 
	krige_data[:, 3] = krige_data[:, 3] .+ sample_f_norm_mean;  # we note that variance is translation invariant

	# make all land 0 (as is required by the model) ::
	krige_data[findall(==(0), query_mask[:]), 3] .= 0; 
	krige_data[findall(==(0), query_mask[:]), 4] .= 0; 

	# reshape and plot ::
	krige_data_grid = reshape(krige_data[:, 3], size(query_mask_grid)); 
	krige_var_grid = reshape(krige_data[:, 4], size(query_mask_grid));

	# save data ::
	pmip3_data[:, :, imonth] = krige_data_grid; 		
	pmip3_var[:, :, imonth] = krige_var_grid; 		

	# plot outputs ::
	fig = plt.figure(figsize=[14, 8]);
	gs = fig.add_gridspec(2, 4, wspace=1.0); 
	ax1 = fig.add_subplot(py"$(gs)[0, 1:3]");
	ax2 = fig.add_subplot(py"$(gs)[1, :2]");
	ax3 = fig.add_subplot(py"$(gs)[1, 2:]");
	cmap = plt.cm.cividis;
	norm = colors.BoundaryNorm(0. : 0.1 : 1.0, cmap.N);
	norm_var = colors.BoundaryNorm(0. : 0.001 : 0.01, cmap.N);
	im1 = ax1.scatter(sample_x_grid, sample_y_grid, c=sample_f, s=50, cmap=cmap, norm=norm);
	divider1 = axes_grid1.make_axes_locatable(ax1);
	cax1 = divider1.append_axes("right", size="5%", pad=0.05);
	fig.colorbar(im1, cax=cax1, orientation="vertical", cmap="cividis");
	ax1.set_title("Sample Data");
	ax1.set_xlim([0, 360]); 
	ax1.set_ylim([-90, 90]); 
	# im2 = ax2.scatter(krige_data[:, 1], krige_data[:, 2], c=krige_data[:, 3], s=150, cmap=cmap, norm=norm);
	im2 = ax2.contourf(query_x_grid, query_y_grid, krige_data_grid, cmap=cmap, norm=norm);
	divider2 = axes_grid1.make_axes_locatable(ax2);
	cax2 = divider2.append_axes("right", size="5%", pad=0.05);
	fig.colorbar(im2, cax=cax2, orientation="vertical", cmap="cividis");
	ax2.set_title("Kriged Data"); 
	ax2.set_xlim([0, 360]); 
	ax2.set_ylim([-90, 90]); 
	# im3 = ax3.scatter(krige_data[:, 1], krige_data[:, 2], c=krige_data[:, 4], s=150, cmap=cmap, norm=norm_var);
	im3 = ax3.contourf(query_x_grid, query_y_grid, krige_var_grid, cmap=cmap, norm=norm_var);
	divider3 = axes_grid1.make_axes_locatable(ax3);
	cax3 = divider3.append_axes("right", size="5%", pad=0.05);
	fig.colorbar(im3, cax=cax3, orientation="vertical", cmap="cividis");
	ax3.set_title("Kriged Variance"); 
	ax3.set_xlim([0, 360]); 
	ax3.set_ylim([-90, 90]);
	plt.savefig(joinpath(basepath_plots, model_name, model_name * "_" * string(imonth) * "_krige_data.png"), dpi=300); 
	plt.show();

	# print out ::
	println("Done kriging data.");

	# close all plots ::
	plt.close("all");

end

# export data to .nc file ::
pmip3_out = NCDataset(joinpath(basepath_krige, model_name, model_name * "_krige.nc"), "c"); 
defDim(pmip3_out, "lon", size(pmip3_data)[1]);
defDim(pmip3_out, "lat", size(pmip3_data)[2]);
defDim(pmip3_out, "time", size(pmip3_data)[3]);
pmip3_out.attrib["title"] = model_name * " SIC Kriged to UVic Grid"; 
sic_data = defVar(pmip3_out, "sic", Float64, ("lon", "lat", "time"));
sic_var = defVar(pmip3_out, "sic_var", Float64, ("lon", "lat", "time"));
sic_data[:, :, :] = pmip3_data;
sic_var[:, :, :] = pmip3_var;
sic_data.attrib["units"] = "sea ice fraction [%]";
sic_var.attrib["units"] = "sea ice fraction variance [%^2]";
close(pmip3_out);

# =========================================================
