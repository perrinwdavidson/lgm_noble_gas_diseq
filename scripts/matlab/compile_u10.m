%%=========================================================================
%   compile_u10
%%-------------------------------------------------------------------------
%   purpose: to compile cmip5 meridional wind speed (u10).
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 01.07.22
%%=========================================================================
%%  configure
%   set the important paths that we will use ::
cmip5_path = '/Volumes/data/CMIP5/output1/';

%   set variable paths ::
lgm_variable_path = 'lgm/mon/atmos/Amon/r1i1p1/ua/';
lgm_paths = {cmip5_path, lgm_variable_path};
pic_variable_path = 'piControl/mon/atmos/Amon/r1i1p1/ua/';
pic_paths = {cmip5_path, pic_variable_path};

%   choose variables ::
variable = 'ua';

%%  read data
cmip5_u10_lgm_raw_data_monthly = read_cmip5_variable(variable, products, lgm_paths);
cmip5_u10_pic_raw_data_monthly = read_cmip5_variable(variable, products, pic_paths); 

%%  make netcdf
write_cmip5_variable(cmip5_u10_lgm_raw_data_monthly, variable, products, 'lgm', fullfile(exp_raw_path, 'cmip5', 'lgm', 'cmip5_u10_lgm_raw_data_monthly.nc'));
write_cmip5_variable(cmip5_u10_pic_raw_data_monthly, variable, products, 'pic', fullfile(exp_raw_path, 'cmip5', 'pic', 'cmip5_u10_pic_raw_data_monthly.nc'));

%%  end program
%%=========================================================================