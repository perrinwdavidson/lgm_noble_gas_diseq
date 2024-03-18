%%=========================================================================
%   compile_masks
%%-------------------------------------------------------------------------
%   purpose: to compile cmip5 land area fraction (sftlf).
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 01.07.22
%%=========================================================================
%%  configure
%   set the important paths that we will use ::
cmip5_path = '/Volumes/data/CMIP5/output1/';

%   set variable paths ::
lgm_variable_path = 'lgm/fx/atmos/fx/r0i0p0/sftlf';
lgm_paths = {cmip5_path, lgm_variable_path};
pic_variable_path = 'piControl/fx/atmos/fx/r0i0p0/sftlf';
pic_paths = {cmip5_path, pic_variable_path};

%   choose variable ::
variable = 'sftlf';

%%  read data
cmip5_sic_lgm_raw_data_monthly = read_cmip5_variable(variable, products, lgm_paths);
cmip5_sic_pic_raw_data_monthly = read_cmip5_variable(variable, products, pic_paths); 

%%  make netcdf
write_cmip5_variable(cmip5_sic_lgm_raw_data_monthly, variable, products, 'lgm', fullfile(exp_raw_path, 'cmip5', 'lgm', 'cmip5_sftlf_lgm_raw_data.nc'))
write_cmip5_variable(cmip5_sic_pic_raw_data_monthly, variable, products, 'pic', fullfile(exp_raw_path, 'cmip5', 'pic', 'cmip5_sftlf_pic_raw_data.nc'))

%%  end program
%%=========================================================================