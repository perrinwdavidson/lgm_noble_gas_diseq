%% ========================================================
%  export_sic
%  --------------------------------------------------------
%  purpose :: to qc and export sic data for interpolation
%  author :: perrin w. davidson
%  contact :: perrinwdavidson@gmail.com
%  data ::26.01.23
%% ========================================================
%% configure
close all;
clear;
clc;

%% load data
%  set basepaths ::
input_path = "data/exp_raw/cmip5/pic/";
output_path = "data/exp_pro/cmip5/sic/pic/";

%  cmip5 data ::
fname = [input_path, "cmip5_sic_pic_raw_data_monthly.nc"];
gnames = ncread(fname, "group_names");
vnames = ncread(fname, "variable_names");

%  uvic data ::
load(fullfile("data", "exp_raw", "uvic", "grid"), "x", "y", "ideep");

%% make variables 
%  uvic mask and grids ::
ideep(ideep ~= 0) = 1; 
grid_mask = logical(ideep);
[grid_lat, grid_lon] = meshgrid(y, x);

%  number of models ::
NUMMODS = length(gnames) + 1;
NUMMONS = 12;

%% initialize netcdf file
fname_export = [output_path, "cmip5_preinterp_data.nc"];

%% loop through all models and export
for iData = 1 : 1 : NUMMODS

	%  only for pmip3 ::
	if iData < 8

		%  load data ::
	    	model_data = ncread(fname, append(gnames{iData}, vnames{iData}));
	        model_lon = ncread(fname, append(gnames{iData}, "lon"));
    		model_lat = ncread(fname, append(gnames{iData}, "lat"));
    
		%  correct values to be in proper domain [0, 1] if needed ::
    		if mean(model_data, "all", "omitnan", "true") > 1
	       		model_data = mean(model_data, 3) ./ 100;
    		end
    
		%    pre-allocate array ::
		model_data_mm = NaN(size(model_data, 1), size(model_data, 2), NUMMON);
        	
		%    loop through all months ::
		for iMonth = 1 : 1 : NUMMON

            		%   calculate mean ::
            		model_data_mm(:, :, iMonth) = mean(model_data(:, :, iMonth : 12 : end), 3, 'omitnan');

        	end
		
		%  make arrays ::
    		model = [model_lon(:), model_lat(:), model_data_mm(:)];
    
		%  remove nans ::
		model_new = model(~isnan(model(:, 3)), :);
    
		%  average duplicate data ::
    		[uniq_coords, ~, idx] = unique(model_new(:, 1:2), "rows");
    		uniq_val = accumarray(idx, model_new(:, 3), [], @mean); 
   	 	model_final = [uniq_coords, uniq_val];

		%  put in netcdf file ::


	%  for uvic data ::
	elseif idata == 8


	end

end

%% export netcdf file 




%% ========================================================
