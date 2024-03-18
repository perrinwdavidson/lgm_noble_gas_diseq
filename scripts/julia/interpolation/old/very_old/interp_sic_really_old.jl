# =========================================================
# interp_sic
# ---------------------------------------------------------
# purpose :: to interpolate sic data to the uvic grid
# author :: perrin w. davidson
# date ::17.02.23
# contact :: perrinwdavidson@gmail.com
# =========================================================
# configure -----------------------------------------------
# load geokrige ::
include("/Users/perrindavidson/.julia/dev/GeoKrige/src/GeoKrige.jl");

# load packages ::
using .GeoKrige  # for kriging
using NCDatasets  # for netcdf functionality
using DelimitedFiles  # for loading delimited files
using Statistics  # for stat functions
using Random  # for random functions
using PyPlot  # for plotting
using PyCall  # for calling python

# specify plotting ::
plt.rc("text", usetex=true);
plt.rc("font", family="serif", size=16);
axes_grid1 = pyimport("mpl_toolkits.axes_grid1");
plt_colors = pyimport("matplotlib.colors"); 

# load data and make inputs -------------------------------
# sic data ::
sic_fname = "data/exp_raw/cmip5/pic/cmip5_sic_pic_raw_data_monthly.nc"; 
sic_nc = NCDataset(sic_fname, "r");
sic_group_names = sic_nc["group_names"][:];
sic_var_names = sic_nc["variable_names"][:];

# query grid ::
query_fname = "data/exp_raw/uvic/grid/";
query_x_grid, query_y_grid, query_mask_grid = readdlm(query_fname * "grid_lon.csv", ',', Float64), 
					      readdlm(query_fname * "grid_lat.csv", ',', Float64),
					      BitArray(readdlm(query_fname * "grid_mask.csv", ',', Float64));

# flatten query grid ::
query_x = query_x_grid[:]; 
query_y = query_y_grid[:]; 
query_mask = query_mask_grid[:]; 
query = hcat(query_x, query_y); 

# set bounds for ice ::
ice_lat_bound = -60;

# dimension of arrays ::
num_models = length(sic_group_names); 
num_month = 12; 

# choose whether or not we will use cressie estimator ::
cressie = false;

# choose whether or not to randomly sample for variogram ::
random_sample = true; 
rand_samp = 0.25; 

# loop through all models ---------------------------------
# loop through models ::
for imod in 1 : 1 : num_models

	# display ::
	println("==== Interpolating " * strip(sic_group_names[imod], '/') * " ===="); 

	# get dataframes ::
	sic_model = sic_nc.group[strip(sic_group_names[imod], '/')];

	# get sic data ::
	sic_data_grid_full = Array{Float64, 3}(sic_model[strip(sic_var_names[imod], '/')][:, :, :]);
	sic_lon_grid = Matrix{Float64}(sic_model["lon"][:, :]);
	sic_lat_grid = Matrix{Float64}(sic_model["lat"][:, :]);

	# make climatology ::
	sic_data_grid = fill!(zeros(size(sic_data_grid_full, 1), size(sic_data_grid_full, 2), num_month), NaN); 
	for imon in 1 : 1 : num_month

		# average data ::
		sic_data_grid[:, :, imon] = mean(sic_data_grid_full[:, :, imon : (num_month + 1) : end], dims=3);

	end

	# if average if greater than 1, ensure that in [0, 1] range ::
	mean_sic = mean(sic_data_grid[findall(!isnan, sic_data_grid)]); 
	if mean_sic > 1
		sic_data_grid /= 100;
	end

	# make plotting directory if does not exist ::
	directory_name = replace(strip(sic_group_names[imod], '/'), " " => "_");
	if !isdir("plots/sic/pic/" * directory_name)
		mkdir("plots/sic/pic/" * directory_name);
	end
	outname = "plots/sic/pic/" * directory_name;

	# make final array ::
	mod_krige_data = fill!(zeros(size(query_mask_grid)[1], size(query_mask_grid)[2], 0), NaN);
	mod_krige_var = fill!(zeros(size(query_mask_grid)[1], size(query_mask_grid)[2], 0), NaN);

	# loop through months ::
	for imon in 1 : 1 : num_month
	
		# artificially impose ice --------------------------------
		# get data ::
		sample_f_grid = sic_data_grid[:, :, imon];

		# find nans ::
		idx_nan = findall(isnan, sample_f_grid);

		# find above lat degrees ::
		idx_lat_bound = findall(<=(ice_lat_bound), sic_lat_grid);  

		# find intersection of both ::
		idx_ice = intersect(idx_lat_bound, idx_nan); 

		# impose ice ::
		sample_f_grid[idx_ice] .= 1.;

		# plot --------------------------------------------------
		plt.scatter(sic_lon_grid, sic_lat_grid, c=sample_f_grid); 
		plt.colorbar(); 
		plt.show();

		# make data ---------------------------------------------
		sample_x, sample_y, sample_f = sic_lon_grid[:], sic_lat_grid[:], sample_f_grid[:];
		
		# quality control data ------------------------------------
		# possibly make longitude [0, 360] ::
		sample_x = shift_lon(sample_x, sample_x);
		
		# get nan idx ::
		idx_sample_real = findall(!isnan, sample_f); 
		
		# select only random subset for data ::
		if random_sample
			idx_sample_real = randsubseq(idx_sample_real, rand_samp)
		end

		# make sample array ::
		sample = hcat(sample_x[idx_sample_real], sample_y[idx_sample_real], sample_f[idx_sample_real]);

		# remove duplicate data :: <- I AM HERE.

		# make inputs ---------------------------------------------
		# normalize and detrend data ::
		sample[:, 3], sample_f_norm_mean = detrend(sample[:, 3]); 

		# plot normalized data ::
		fig, ax = plt.subplots(1, 2, figsize=[12, 6]);
		ax[1].hist(sample_f, bins=sturges(length(sample_f))); 
		ax[1].set_title("Raw");
		ax[1].set_ylabel("Counts");
		ax[1].set_xlabel("Data");
		ax[2].hist(sample[:, 3], bins=sturges(length(sample[:, 3]))); 
		ax[2].set_title("Detrended (Constant)");
		ax[2].set_ylabel("Counts");
		ax[2].set_xlabel("Data");
		plt.show(); 

		# print out ::
		println("Done making inputs.");

		# calculate experimental variogram ------------------------
		# set inputs ::
		max_lag = 5E6; 
		lags, semivar = Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}();

		# calculate matheron estimator ::
		ilags, isemivar = semivariogram(sample[:, 1:2], 
						sample[:, 3], 
						max_lag; 
						estimator="matheron");
		lags = push!(lags, "matheron" => ilags); 
		semivar = push!(semivar, "matheron" => isemivar); 

		# calculate cressie estimator ::
		if cressie
			ilags, isemivar = semivariogram(sample[:, 1:2], 
							sample[:, 3], 
							max_lag; 
							estimator="cressie");
			lags = push!(lags, "cressie" => ilags); 
			semivar = push!(semivar, "cressie" => isemivar); 
		end

		# plot ::
		fig, ax = plt.subplots(1, 1, figsize=[6, 6]);
		ax.scatter(lags["matheron"], semivar["matheron"], label="Matheron");
		if cressie
			ax.scatter(lags["cressie"], semivar["cressie"], label="Cressie");
		end
		ax.set_xlabel("Lag [m]");
		ax.set_ylabel("Semivariance"); 
		ax.set_title("Semivariogram");
		plt.legend(); 
		plt.show();

		# fit semivariogram ---------------------------------------
		# make inputs ::
		vario_estimator = "matheron"; 
		p0 = [0., 1E-1, 2.5E6];

		# choose models ::
		models = ["Spherical", "Gaussian", "Exponential", "Power"]; 

		# fit variogram :: 
		best_model, best_model_name, vario_err, pf, pferr, chisq, dof = fit_variogram(lags[vario_estimator], 
											      semivar[vario_estimator], 
											      p0, 
											      models; 
											      plot_output=true);		

		# display ::
		println("Best model is: ", best_model_name); 

		# krige ---------------------------------------------------
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

		# reshape ::
		krige_data_grid = reshape(krige_data[:, 3], size(query_mask_grid)); 
		krige_var_grid = reshape(krige_data[:, 4], size(query_mask_grid));

		# store ::
		mod_krige_data = cat(mod_krige_data, krige_data_grid, dims=3); 
		mod_krige_var = cat(mod_krige_var, krige_var_grid, dims=3); 

		# plot outputs ::
		fig = plt.figure(figsize=[14, 8]);
		gs = fig.add_gridspec(2, 4, wspace=1.0); 
		ax1 = fig.add_subplot(py"$(gs)[0, 1:3]");
		ax2 = fig.add_subplot(py"$(gs)[1, :2]");
		ax3 = fig.add_subplot(py"$(gs)[1, 2:]");
		cmap = plt.cm.cividis;
		norm = plt_colors.BoundaryNorm(0. : 0.1 : 1.0, cmap.N);
		norm_var = plt_colors.BoundaryNorm(0. : 0.001 : 0.01, cmap.N);
		im1 = ax1.scatter(sample_x, sample_y, c=sample_f, s=50, cmap=cmap, norm=norm);
		divider1 = axes_grid1.make_axes_locatable(ax1);
		cax1 = divider1.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im1, cax=cax1, orientation="vertical", cmap="cividis");
		ax1.set_title("Sample Data");
		ax1.set_xlim([0, 360]); 
		ax1.set_ylim([-90, 90]); 
		im2 = ax2.scatter(krige_data[:, 1], krige_data[:, 2], c=krige_data[:, 3], s=150, cmap=cmap, norm=norm);
		divider2 = axes_grid1.make_axes_locatable(ax2);
		cax2 = divider2.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im2, cax=cax2, orientation="vertical", cmap="cividis");
		ax2.set_title("Kriged Data"); 
		ax2.set_xlim([0, 360]); 
		ax2.set_ylim([-90, 90]); 
		im3 = ax3.scatter(krige_data[:, 1], krige_data[:, 2], c=krige_data[:, 4], s=150, cmap=cmap, norm=norm_var);
		divider3 = axes_grid1.make_axes_locatable(ax3);
		cax3 = divider3.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im3, cax=cax3, orientation="vertical", cmap="cividis");
		ax3.set_title("Kriged Variance"); 
		ax3.set_xlim([0, 360]); 
		ax3.set_ylim([-90, 90]);
		fig.suptitle(strip(sic_group_names[imod], '/') * ": Month " * string(imon)); 
		plt.savefig(outname * "/" * directory_name * "_Month_" * string(imon) * ".png", dpi=300); 
		plt.show();
		plt.close("all"); 

	end

	# make data directory if does not exist ::
	directory_name = replace(strip(sic_group_names[imod], '/'), " " => "_");
	if !isdir("data/sims/cmip5/sic/pic/" * directory_name)
		mkdir("data/sims/cmip5/sic/pic/" * directory_name);
	end
	outname = "data/sims/cmip5/sic/pic/" * directory_name * "/" * directory_name * "_sic.nc";

	# save data ::
	data_out = NCDataset(outname, "c");
	defDim(data_out, "lon", size(query_x_grid)[1]);
	defDim(data_out, "lat", size(query_y_grid)[2]);
	defDim(data_out, "time", num_month); 
	data_out.attrib["title"] = directory_name * " SIC"; 
	sic_data_out = defVar(data_out, "SIC", Float64, ("lon", "lat", "time"));
	sic_var_out = defVar(data_out, "Variance", Float64, ("lon", "lat", "time"));
	sic_lon = defVar(data_out, "lon", Float64, ("lon", "lat"));
	sic_lat = defVar(data_out, "lat", Float64, ("lon", "lat"));
	sic_data_out[:, :, :] = mod_krige_data[:, :, :];
	sic_var_out[:, :, :] = mod_krige_var[:, :, :];
	sic_lon[:, :] = query_x_grid; 
	sic_lat[:, :] = query_y_grid; 
	sic_data_out.attrib["units"] = "% ice coverage"
	sic_var_out.attrib["units"] = "[% ice coverage]^2"
	close(data_out)
	
end

# =========================================================
