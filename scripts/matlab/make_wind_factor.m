%%=========================================================================
%   make_wind_factor
%%-------------------------------------------------------------------------
%   purpose: to make a pic to lgm wind factor for core-2 winds.
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 28.07.22
%%=========================================================================
%%  get core grid
%   get vectors
x = ncread(fullfile('data', 'exp_raw', 'core2', 'u_10.15JUNE2009.nc'), 'LON');
y = ncread(fullfile('data', 'exp_raw', 'core2', 'u_10.15JUNE2009.nc'), 'LAT');
t = ncread(fullfile('data', 'exp_raw', 'core2', 'u_10.15JUNE2009.nc'), 'TIME');

%   make grid ::
[core_lat, core_lon] = meshgrid(y, x);  

%   get lengths ::
nx = length(x); 
ny = length(y);
nt = length(t);

%%   make initial data
%    set wind factors ::
wind_factors = 0.5 : 0.1 : 1.5; 
NUMFACTOR = length(wind_factors); 

%    set extent ::
lat_extent_min = -50; 
lat_extent_max = 50; 

%    base ones array ::
one_array = ones(size(core_lat));
one_array = repmat(one_array, [1 1 nt]); 

%    pre-allocate final array ::
wind_factor_array = cell(2, NUMFACTOR);  % (1, :) is full extent, (2, :) is HL 
variable_names = cell(1, NUMFACTOR); 

%%   loop through all data
for iFactor = 1 : 1 : NUMFACTOR

	%    get variable name ::
	variable_name = num2str(wind_factors(iFactor)); 
	variable_name = erase(variable_name, '.');
	variable_names{iFactor} = variable_name;

	%    make factor arrays ::
	%%%  full extent ::
	wind_factor_array{1, iFactor} = wind_factors(iFactor) .* one_array; 

	%%%  hl ::
	wind_factor_array{2, iFactor} = wind_factors(iFactor) .* one_array; 
	hl_idx = find( (y > lat_extent_min) & (y < lat_extent_max));
	wind_factor_array{2, iFactor}(:, hl_idx, :) = 1;  

end

%%   write data to netcdf4 files
%    re-write ::
delete(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'));
delete(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'));

%    create ::
%%%  lon ::
nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), 'lon', 'dimensions', {'lon', nx});
nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), 'lon', 'dimensions', {'lon', nx});

%%%  lat ::
nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), 'lat', 'dimensions', {'lat', ny});
nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), 'lat', 'dimensions', {'lat', ny});

%%%  time ::
nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), 'time', 'dimensions', {'time', nt});
nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), 'time', 'dimensions', {'time', nt});

%%%  factor ::
for iFactor = 1 : 1 : NUMFACTOR

	%     create ::
	nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), append('windfactor_', variable_names{iFactor}), 'dimensions', {'lon', nx, 'lat', ny, 'time', nt});
	nccreate(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), append('windfactor_', variable_names{iFactor}), 'dimensions', {'lon', nx, 'lat', ny, 'time', nt});

end

%    write ::
%%%  lon ::
ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), 'lon', x);
ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), 'lon', x);

%%%  lat ::
ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), 'lat', y);
ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), 'lat', y);

%%%  time ::
ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), 'time', t);
ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), 'time', t);

%%%  factor ::
for iFactor = 1 : 1 : NUMFACTOR

	%     write ::
	ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_90S90N.nc'), append('windfactor_', variable_names{iFactor}), wind_factor_array{1, iFactor});
	ncwrite(fullfile('data', 'sims', 'wind_factor', 'wind_factor_50S50N.nc'), append('windfactor_', variable_names{iFactor}), wind_factor_array{2, iFactor});

end

%%=========================================================================
