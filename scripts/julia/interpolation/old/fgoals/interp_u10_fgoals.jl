# =========================================================
# interp_u10_fgoals
# ---------------------------------------------------------
# purpose :: to interpolate fgoals u10 data to uvic grid
# author :: perrin w. davidson
# contact :: perrinwdavidson@gmail.com
# date :: 05.04.23
# =========================================================
function interp_u10_fgoals(fname, basepath_sample, basepath_plots, basepath_interp)

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
	sample_x_raw, sample_y_raw = sample_x_grid[:], sample_y_grid[:];

	# make longitude [0, 360] ::
	sample_x_raw = shift_lon(sample_x_raw, sample_x_raw); 
	
	# make coordinate array ::
	sample_coords = hcat(sample_x_raw, sample_y_raw); 

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

	# set ice bound ::
	ice_lat_bound = -65;

	# make directory if it doesn't exist ::
	if !isdir(joinpath(basepath_plots, model_name))
		mkdir(joinpath(basepath_plots, model_name)); 
	end
	if !isdir(joinpath(basepath_interp, model_name))
		mkdir(joinpath(basepath_interp, model_name)); 
	end

	# preallocate final array ::
	pmip3_data = zeros(size(query_mask_grid)[1], size(query_mask_grid)[2], num_months); 

	# print out ::
	println("Done loading data.");

	# loop through all months ---------------------------------
	for imonth in 1 : 1 : num_months

		# quality control data ----------------------------
		# get monthly data ::
		sample_f_raw_grid = sample_f_grid[:, :, imonth];
		sample_f_raw = sample_f_raw_grid[:];

		# remove duplicate data ::
		idx_duplicate = groupinds(groupslices(sample_coords, dims=1));
		sample_x = [];
		sample_y = [];
		sample_f = [];
		for i in idx_duplicate
			sample_x = vcat(sample_x, sample_coords[i[1], 1]); 
			sample_y = vcat(sample_y, sample_coords[i[1], 2]); 
			sample_f = vcat(sample_f, mean(sample_f_raw[i])); 
		end
		
		# get nan idx ::
		idx_sample_real = findall(!ismissing, sample_f[:]); 
		
		# make input data ::
		sample = Matrix{Float64}(hcat(sample_x[idx_sample_real], sample_y[idx_sample_real], sample_f[idx_sample_real])); 
		points = sample[:, 1:2]';
		samples = sample[:, 3];
		query_points = query';
		
		# print out ::
		println("Done making inputs.");

		# make interpolator ::
		interpolator = ScatteredInterpolation.interpolate(NearestNeighbor(), points, samples);

		# interpolate data ::
		interp_data = ScatteredInterpolation.evaluate(interpolator, query_points); 
		
		# make all land 0 (as is required by the model) ::
		interp_data[findall(==(0), query_mask[:])] .= 0; 

		# reshape and plot ::
		interp_data_grid = reshape(interp_data, size(query_mask_grid)); 
		
		# print out ::
		println("Done interpolating data.");

		# plot outputs ::
		fig, ax = plt.subplots(1, 2, figsize=[18, 6]); 
		cmap = plt.cm.cividis;
		norm = colors.BoundaryNorm(-20: 2 : 20, cmap.N);
		im1 = ax[1].scatter(sample_x[idx_sample_real], sample_y[idx_sample_real], c=sample_f[idx_sample_real], s=50, cmap=cmap, norm=norm);
		divider1 = axes_grid1.make_axes_locatable(ax[1]);
		cax1 = divider1.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im1, cax=cax1, orientation="vertical", cmap="cividis");
		ax[1].set_title("Sample Data");
		ax[1].set_xlim([0, 360]); 
		ax[1].set_ylim([-90, 90]); 
		im2 = ax[2].contourf(query_x_grid, query_y_grid, interp_data_grid, cmap=cmap, norm=norm);
		divider2 = axes_grid1.make_axes_locatable(ax[2]);
		cax2 = divider2.append_axes("right", size="5%", pad=0.05);
		fig.colorbar(im2, cax=cax2, orientation="vertical", cmap="cividis");
		ax[2].set_title("Interpolated Data"); 
		ax[2].set_xlim([0, 360]); 
		ax[2].set_ylim([-90, 90]);
		plt.tight_layout(); 
		plt.savefig(joinpath(basepath_plots, model_name, model_name * "_" * var_name * "_" * string(imonth) * "_interp_data.png"), dpi=300); 
		plt.show();

		# save data ::
		pmip3_data[:, :, imonth] = interp_data_grid; 	

		# print out ::
		println("Done with month: ", imonth);

	end

	# export data to .nc file ::
	pmip3_out = NCDataset(joinpath(basepath_interp, model_name, model_name * "_" * var_name * "_krige.nc"), "c"); 
	defDim(pmip3_out, "lon", size(pmip3_data)[1]);
	defDim(pmip3_out, "lat", size(pmip3_data)[2]);
	defDim(pmip3_out, "time", size(pmip3_data)[3]);
	pmip3_out.attrib["title"] = model_name * " " * uppercase(var_name) * " NN Interpolated to UVic Grid"; 
	wind_data = defVar(pmip3_out, var_name, Float64, ("lon", "lat", "time"));
	wind_data[:, :, :] = pmip3_data;
	wind_data.attrib["units"] = "velocity [m s-1]";
	close(pmip3_out);

	# print out ::
	println("Done with model: ", model_name);

end

# =========================================================
