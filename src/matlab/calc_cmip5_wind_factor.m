function model_wind_factor_mean = calc_cmip5_wind_factor(filename)
%%-------------------------------------------------------------------------
%   purpose: to calculate the cmip5 windspeed weighted mean
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 06.08.22
%%-------------------------------------------------------------------------
%%  configure
%   set two filenames ::
filename_u = filename{1};  % will be the file which we draw everything from, assuming that there is no difference between u and v. 
filename_v = filename{2}; 

%   set bounds ::
lat_min_bound = -50; 
lat_max_bound = 50; 

%   get age ::
age = filename_u(11 : 13);
variable = filename_u(7 : 9);

%   get variables names ::
%%% u ::
group_names_u = ncread(filename_u, 'group_names');
variable_names_u = ncread(filename_u, 'variable_names');

%%% v ::
group_names_v = ncread(filename_v, 'group_names');
variable_names_v = ncread(filename_v, 'variable_names');

%%% general ::
group_names = group_names_u; 
variable_names = variable_names_u; 

%   get land mask filename ::
if strcmp(age, 'lgm')

    mask_filename = 'land_mask_lgm.nc';

elseif strcmp(age, 'pic')

    mask_filename = 'land_mask_pic.nc';

end

%   get land mask variables names ::
group_names_mask = ncread(mask_filename, 'group_names');
variable_names_mask = ncread(mask_filename, 'variable_names');

%%  set up uvic
%   load uvic ::
load(fullfile('data', 'exp_raw', 'uvic', 'grid'), 'x', 'y', 'ideep');
load(fullfile('data', 'exp_raw', 'uvic', age, 'wind_speed'), 'windspeed');
[uvic_lat, uvic_lon] = meshgrid(y, x); 

%    append names ::
group_names = [group_names; {append('UVic ', upper(age), ' Default')}];
variable_names = [variable_names; {append('UVic Default ', upper(age), ' ' , upper(variable))}];

%   make land mask ::
ideep(ideep ~= 0) = 1; 
ind_lm = find(ideep(:) == 1);
uvic_land_mask = logical(ideep(:));

%%  finalize loop variables and functions
%   get number of models ::
NUMMOD = size(group_names, 1);
NUMMON = 12; 

%   set averaging months ::
months_south = [6, 7, 8];  % june, july, august 
months_north = [12, 1, 2];   % december, january, february

%   set days per month ::
days_south = [30, 31, 31]; 
days_north = [31, 31, 28];  % no leap year

%   define wind magnitude function ::
wind_mag = @(u, v) sqrt((u .^ 2) + (v .^ 2));

%%  loop through all data products
for iMod = 1 : 1 : NUMMOD

	if iMod < 8

		%     get data ::
      model_data_u = ncread(filename_u, append(group_names_u{iMod}, variable_names_u{iMod}));
      model_data_v = ncread(filename_v, append(group_names_v{iMod}, variable_names_v{iMod}));
		model_lon = ncread(filename_u, append(group_names_u{iMod}, 'lon'));
      model_lat = ncread(filename_u, append(group_names_u{iMod}, 'lat'));
      model_mask = ncread(mask_filename, append(group_names_mask{iMod}, variable_names_mask{iMod}));

		%     calculate magnitude ::
		model_data = wind_mag(model_data_u, model_data_v); 

	elseif iMod == 8

		%    get data ::
		model_data = windspeed;
		model_lon = uvic_lon; 
		model_lat = uvic_lat;
		model_mask = uvic_land_mask; 

	end

	%   average data ::
    	if iMod < 8 
		
		%   display what is going on ::
        	disp('Calculating monthly mean climatology.');
        
		%    pre-allocate array ::
		model_data_mm = NaN(size(model_data, 1), size(model_data, 2), NUMMON);
        	
		%    loop through all months ::
		for iMonth = 1 : 1 : NUMMON

            		%   calculate mean ::
            		model_data_mm(:, :, iMonth) = mean(model_data(:, :, iMonth : 12 : end), 3, 'omitnan');

        	end

		%    store data ::
        	model_data = model_data_mm;

    	end

	%     get size ::
	[nx, ny, nt] = size(model_data); 
	
	%     mask data ::
	model_data(~repmat(model_mask, [1 1 nt])) = NaN; 	

	%     start sums ::
	%%%   southern hemisphere ::
	data_sum_south = 0; 
	data_weight_south = 0; 
	iMonS = 1; 

	%%%   northern hemisphere ::
	data_sum_north = 0; 
	data_weight_north = 0; 
	iMonN = 1; 

	%     loop through all time ::
	for it = 1 : 1 : nt

		%     see if should count ::
		count_south = 0; 
		count_north = 0; 
		if ismember(it, months_south)

			count_south = 1; 

		elseif ismember(it, months_north)

			count_north = 1; 

		end

		%     loop through all long ::
		for ix = 1 : 1 : nx

			%     loop through all lat ::
			for iy = 1 : 1 : ny

				%     southern hemisphere months ::
				if (~isnan(model_data(ix, iy, it))) && (model_lat(ix, iy) <= lat_min_bound) && ismember(it, months_south)

					%     sum data ::
					data_sum_south = data_sum_south + (model_data(ix, iy, it) .* cosd(model_lat(ix, iy))); 
					data_weight_south = data_weight_south + cosd(model_lat(ix, iy)); 
				
				end

				%     northern hemisphere months ::
				if (~isnan(model_data(ix, iy, it))) && (model_lat(ix, iy) >= lat_max_bound) && ismember(it, months_north)

					%     sum data ::
					data_sum_north = data_sum_north + (model_data(ix, iy, it) .* cosd(model_lat(ix, iy))); 
					data_weight_north = data_weight_north + cosd(model_lat(ix, iy)); 

				end

			end

		end

		%    iterate and calculate weighted mean ::
		if count_south

			model_wind_factor.south(iMonS, iMod) = data_sum_south / data_weight_south; 
			iMonS = iMonS + 1; 
		
		elseif count_north

			model_wind_factor.north(iMonN, iMod) = data_sum_north / data_weight_north; 
			iMonN = iMonN + 1; 
		
		end

	end

	%    display what is going on ::
	disp(append('Done with ', upper(age), ' model ', num2str(iMod))); 

end

%%  calculate seasonal day weighted average 
%   make weight arrays ::
day_wind_factor.north = repmat(days_north', [1 NUMMOD]) .* model_wind_factor.north; 
day_wind_factor.south = repmat(days_south', [1 NUMMOD]) .* model_wind_factor.south; 
total_days.north = sum(days_north, 'all'); 
total_days.south = sum(days_south, 'all'); 

%   calculate seasonal weighted average ::
model_wind_factor_mean.north = sum(day_wind_factor.north) ./ total_days.north; 
model_wind_factor_mean.south = sum(day_wind_factor.south) ./ total_days.south; 

%%-------------------------------------------------------------------------
end
