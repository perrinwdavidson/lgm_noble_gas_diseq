%%=========================================================================
%   run_cmip5
%%-------------------------------------------------------------------------
%   purpose: to compile cmip5 windspeed (u10, v10) and sea ice fraction (sic).
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 01.07.22
%%=========================================================================
%%  configure
%   deal with the environment ::
clc;
cd('/Users/perrindavidson/research/whoi/current/noble_gas/cmip5')
clear; 
close;

%   add all subdirectories to path ::
addpath(genpath(pwd));

%   set input and output paths ::
exp_raw_path = fullfile(pwd, 'data', 'exp_raw'); 
exp_pro_path = fullfile(pwd, 'data', 'exp_pro');

%%  set global model parameters
%   write down the products that you want of the form {institute, model} ::
products = {'NCAR', 'CCSM4'; ...
            'MRI', 'MRI-CGCM3'; ...
            'LASG-CESS', 'FGOALS-g2'; ...
            'IPSL', 'IPSL-CM5A-LR'; ...
            'MPI-M', 'MPI-ESM-P'; ...
            'CNRM-CERFACS', 'CNRM-CM5'; ...
            'MIROC', 'MIROC-ESM'};

%%  compile data
%   sea ice fraction ::
compile_sic; 

%   windspeed ::
compile_u10; 
compile_v10;

%   land masks ::
compile_masks;  

%%  make land mask
make_masks; 

%%  make wind factor;
make_wind_factor;

%%  calculate cmip5 ratios
calc_cmip5_wind_factor;  % <- I am here.

%%  plot output
plot_tmm_output;

%%  end program
%%=========================================================================
