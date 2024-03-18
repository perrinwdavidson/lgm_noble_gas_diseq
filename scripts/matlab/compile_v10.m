%%=========================================================================
%   compile_v10
%%-------------------------------------------------------------------------
%   purpose: to compile cmip5 zonal wind speed (v10).
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 01.07.22
%%=========================================================================
%%  configure
%   set the important paths that we will use ::
cmip5_path = '/Volumes/data/CMIP5/output1/';

%   set variable paths ::
lgm_variable_path = 'lgm/mon/atmos/Amon/r1i1p1/va/';
lgm_paths = {cmip5_path, lgm_variable_path};
pic_variable_path = 'piControl/mon/atmos/Amon/r1i1p1/va/';
pic_paths = {cmip5_path, pic_variable_path};

%   choose variables ::
variable = 'va';

%%  read data
cmip5_v10_lgm_raw_data_monthly = read_cmip5_variable(variable, products, lgm_paths);
cmip5_v10_pic_raw_data_monthly = read_cmip5_variable(variable, products, pic_paths); 

%%  make netcdf
write_cmip5_variable(cmip5_v10_lgm_raw_data_monthly, variable, products, 'lgm', fullfile(exp_raw_path, "cmip5", "lgm", "cmip5_v10_lgm_raw_data_monthly.nc"));
write_cmip5_variable(cmip5_v10_pic_raw_data_monthly, variable, products, 'pic', fullfile(exp_raw_path, "cmip5", "pic", "cmip5_v10_pic_raw_data_monthly.nc"));

%%  end program
%%=========================================================================