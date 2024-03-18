%%=========================================================================
%   estimate_wind_distribution
%%-------------------------------------------------------------------------
%   purpose: to estimate the cmip5 wind distribution for high latitudes
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 22.09.22
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
cmip5_wind_distribution = cell(1, length(filenames)); 

%   loop through all files ::
for iFile = filenames

    %   get filename ::
    filename = iFile{:};

    %   calculate factor the file ::
    cmip5_wind_distribution{iCount} = calc_cmip5_wind_distribution(filename);

    %   count ::
    iCount = iCount + 1; 

end
