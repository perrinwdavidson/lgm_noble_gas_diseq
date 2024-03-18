%%=========================================================================
%   make_cmip5_wind_factor
%%-------------------------------------------------------------------------
%   purpose: to calculate the pmip3 wind factor for high latitudes
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 06.08.22
%%=========================================================================
%%  configure
%   set file names ::
filenames = {{'cmip5_u10_lgm_raw_data_monthly.nc', ...
             'cmip5_v10_lgm_raw_data_monthly.nc'}, ...
             {'cmip5_u10_pic_raw_data_monthly.nc', ...
             'cmip5_v10_pic_raw_data_monthly.nc'}};  % note {u, v} and {lgm, pic} filenaming convention 

%%  calculate individual ages
%   start counter ::
iCount = 1; 

%   pre-allocate ::
cmip5_wind_factor = cell(1, length(filenames)); 

%   loop through all files ::
for iFile = filenames

    %   get filename ::
    filename = iFile{:};

    %   calculate factor the file ::
    cmip5_wind_factor{iCount} = calc_cmip5_wind_factor(filename);

    %   count ::
    iCount = iCount + 1; 

end

%%  calculate final factors
%   calculate ratio ::
wind_factor_north = cmip5_wind_factor{1}.north ./ cmip5_wind_factor{2}.north;
wind_factor_south = cmip5_wind_factor{1}.south ./ cmip5_wind_factor{2}.south;

%%  save data
%   make array ::
wind_factor.north = wind_factor_north; 
wind_factor.south = wind_factor_south; 

%   save array ::
save(fullfile('data', 'sims', 'cmip5', 'wind_factor', 'wind_factor'), 'wind_factor'); 

%%   end program
