# =========================================================
# interp_u10_ccsm4
# ---------------------------------------------------------
# purpose :: to interpolate ccsm4 u10 data to uvic grid
# author :: perrin w. davidson
# contact :: perrinwdavidson@gmail.com
# date :: 05.04.23
# =========================================================
function interp_u10_ccsm4(fname, basepath_sample, basepath_plots, basepath_krige)

	# load data -----------------------------------------------
	# get names ::
	model_name = lowercase(split(split(fname, "_")[3], "-")[1]);
	var_name = split(fname, "_")[1];

	# sample data ::
	sample_nc = NCDataset(joinpath(basepath_sample, fname), "r");
	sample_x_grid = Array{Float64}(sample_nc["lon"][:]);
	sample_y_grid = Array{Float64}(sample_nc["lat"][:]);
	sample_f_grid = sample_nc[var_name][:, :, 1, :];
	if length(size(sample_x_grid)) == 1
		sample_x_grid, sample_y_grid = meshgrid(sample_x_grid, sample_y_grid);
	end
	sample_x, sample_y = sample_x_grid[:], sample_y_grid[:];

	# make longitude [0, 360] ::
	sample_x = shift_lon(sample_x, sample_x); 

	# query grid ::
	query_x_grid = readdlm(joinpath(basepath_query, "grid_lon.csv"), ',', Float64);
	query_y_grid = readdlm(joinpath(basepath_query, "grid_lat.csv"), ',', Float64);
	query_mask_grid = BitArray(readdlm(joinpath(basepath_query, "grid_mask.csv"), ',', Float64)); 
		
	# flatten query grid ::
	query_x = query_x_grid[:]; 
	query_y = query_y_grid[:]; 
	query_mask = query_mask_grid[:]; 
	query = hcat(query_x, query_y); 

	# get number of months ::
	num_months = size(sample_f_grid)[3];

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

		# get nan idx ::
		idx_sample_real = findall(!ismissing, sample_f[:]); 

		# make inputs -------------------------------------
		# set normalization and detrending ::
		norm_data = false;
		detrend_data = false;

		# get constant to add to non-negative values and normalize data ::
		if norm_data
			pos_constant = ceil(abs(minimum(sample_f[idx_sample_real]))); 
			sample_f_norm, lambda_f = normalize(sample_f[idx_sample_real] .+ pos_constant); 
		else
			sample_f_norm = Array{Float64}(sample_f[idx_sample_real]);
		end

		# detrend data ::
		if detrend_data
			sample_f_norm_detrend, sample_f_norm_mean = detrend(sample_f_norm); 
		else
			sample_f_norm_detrend = sample_f_norm; 
		end

		# make sample array ::
		sample = hcat(sample_x[idx_sample_real], sample_y[idx_sample_real], sample_f_norm_detrend); 

		# plot normalized data ::
		fig, ax = plt.subplots(1, 2, figsize=[12, 6]);
		ax[1].hist(sample_f[idx_sample_real], bins=sturges(length(sample_f))); 
		ax[1].set_title("Raw");
		ax[1].set_ylabel("Counts");
		ax[1].set_xlabel("Data");
		ax[2].hist(sample_f_norm_detrend, bins=sturges(length(sample_f_norm_detrend))); 
		ax[2].set_title("Normalized");
		ax[2].set_ylabel("Counts");
		ax[2].set_xlabel("Data");
		plt.tight_layout(); 
		plt.savefig(joinpath(basepath_plots, model_name, model_name * "_" * var_name * "_" * string(imonth) * "_norm_detrend_data.png"), dpi=300); 
		plt.show(); 

		# print out ::
		println("Done making inputs.");

		# calculate experimental variogram ----------------
		# set inputs ::
		max_lag = 5E6; 
		lags, semivar = Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}();
		vario_estimator = "cressie"; 
		subset_perc = 0.1;
		plot_exp_vario = false;

		# choose subset of data ::
		idx_subset = randsubseq(MersenneTwister(1234), 1 : 1 : length(sample_f_norm_detrend), subset_perc)

		if vario_estimator == "matheron"

			# calculate matheron estimator ::
			ilags, isemivar = semivariogram(sample[idx_subset, 1:2], 
							sample[idx_subset, 3], 
							max_lag; 
							estimator="matheron");
			lags = push!(lags, "matheron" => ilags); 
			semivar = push!(semivar, "matheron" => isemivar); 

		elseif vario_estimator == "cressie"

			# calculate cressie estimator ::
			ilags, isemivar = semivariogram(sample[idx_subset, 1:2], 
							sample[idx_subset, 3], 
							max_lag; 
							estimator="cressie");
			lags = push!(lags, "cressie" => ilags); 
			semivar = push!(semivar, "cressie" => isemivar); 

		end

		# plot ::
		if plot_exp_vario
		 	fig, ax = plt.subplots(1, 1, figsize=[6, 6]);
		 	ax.scatter(lags["matheron"], semivar["matheron"], label="Matheron");
		 	ax.scatter(lags["cressie"], semivar["cressie"], label="Cressie");
		 	ax.set_xlabel("Lag [m]");
		 	ax.set_ylabel("Semivariance"); 
		 	ax.set_title("Semivariogram");
		 	plt.legend();
		 	plt.tight_layout(); 
			plt.show();
		end

		# fit semivariogram -------------------------------
		# make inputs ::
		p0 = [0., 1E-2, 2E6];

		# choose models ::
		models = ["Spherical"];

		# fit variogram :: 
		best_model, best_model_name, vario_err, pf, pferr, chisq, dof = fit_variogram(lags[vario_estimator], 
											      semivar[vario_estimator], 
											      p0, 
											      models; 
											      plot_output=true,
											      save_fig=joinpath(basepath_plots, model_name, model_name * "_" * var_name * "_" * string(imonth) * "_fit_semivariogram.png"));		

		# display ::
		println("Done calculating and fitting experimental variogram.");
		println("Best model is: ", best_model_name); 

		# krige -------------------------------------------
		# set parameters ::
		max_lag_krige = 1E5; 

		# krige data ::
		krige_data_norm_detrend = krige(sample,
						query,
						query_mask,
						best_model,
						pf[best_model_name],
						max_lag_krige); 	

		# renormalize ::
		if norm_data
			krige_data_norm = deepcopy(krige_data_norm_detrend); 
			krige_data_norm[:, 3] = krige_data_norm[:, 3] .+ sample_f_norm_mean;  # we note that variance is translation invariant
		else
			krige_data_norm = deepcopy(krige_data_norm_detrend);
		end

		# trend data ::
		if detrend_data
			kriged_data = deepcopy(krige_data_norm); 
			for i in [3, 4]
				kriged_data[:, i] = denormalize(krige_data_norm[:, i], lambda_f); 
			end
			kriged_data[:, 3] .=- pos_constant;  # we again note that variance is translation invariance
		else
			kriged_data = krige_data_norm; 
		end

		# make all land 0 (as is required by the model) ::
		kriged_data[findall(==(0), query_mask[:]), 3] .= 0; 
		kriged_data[findall(==(0), query_mask[:]), 4] .= 0; 

		# reshape and plot ::
		krige_data_grid = reshape(kriged_data[:, 3], size(query_mask_grid)); 
		krige_var_grid = reshape(kriged_data[:, 4], size(query_mask_grid));

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
		im1 = ax1.scatter(sample_x[idx_sample_real], sample_y[idx_sample_real], c=sample_f[idx_sample_real], s=50);
		divider1 = axes_grid1.make_axes_locatable(ax1);
		cax1 = divider1.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im1, cax=cax1, orientation="vertical", cmap="cividis");
		ax1.set_title("Sample Data");
		ax1.set_xlim([0, 360]); 
		ax1.set_ylim([-90, 90]); 
		im2 = ax2.contourf(query_x_grid, query_y_grid, krige_data_grid); 
		divider2 = axes_grid1.make_axes_locatable(ax2);
		cax2 = divider2.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im2, cax=cax2, orientation="vertical", cmap="cividis");
		ax2.set_title("Kriged Data"); 
		ax2.set_xlim([0, 360]); 
		ax2.set_ylim([-90, 90]); 
		im3 = ax3.contourf(query_x_grid, query_y_grid, krige_var_grid);
		divider3 = axes_grid1.make_axes_locatable(ax3);
		cax3 = divider3.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im3, cax=cax3, orientation="vertical", cmap="cividis");
		ax3.set_title("Kriged Variance"); 
		ax3.set_xlim([0, 360]); 
		ax3.set_ylim([-90, 90]);
		plt.savefig(joinpath(basepath_plots, model_name, model_name * "_" * var_name * "_" * string(imonth) * "_krige_data.png"), dpi=300); 
		plt.show();

		# print out ::
		println("Done kriging data.");

		# close all plots ::
		plt.close("all");

	end

	# export data to .nc file ::
	pmip3_out = NCDataset(joinpath(basepath_krige, model_name, model_name * "_" * var_name * "_krige.nc"), "c"); 
	defDim(pmip3_out, "lon", size(pmip3_data)[1]);
	defDim(pmip3_out, "lat", size(pmip3_data)[2]);
	defDim(pmip3_out, "time", size(pmip3_data)[3]);
	pmip3_out.attrib["title"] = model_name * " " * uppercase(var_name) * " Kriged to UVic Grid"; 
	wind_data = defVar(pmip3_out, var_name, Float64, ("lon", "lat", "time"));
	wind_var = defVar(pmip3_out, var_name * "_var", Float64, ("lon", "lat", "time"));
	wind_data[:, :, :] = pmip3_data;
	wind_var[:, :, :] = pmip3_var;
	wind_data.attrib["units"] = "speed [m s-1]";
	wind_var.attrib["units"] = "speed variance [m+2 s-2]";
	close(pmip3_out);

end

# =========================================================
