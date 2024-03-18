%%=========================================================================
%   make_masks
%%-------------------------------------------------------------------------
%   purpose: to make masks for cmpi5 and uvic.
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 22.07.22
%%=========================================================================
%%  configure
%   set file names ::
filenames = {'cmip5_sftlf_lgm_raw_data.nc', ...
             'cmip5_sftlf_pic_raw_data.nc'}; 

%   set critical land fraction value below which is ocean ::
land_crit = 0.50; 

%%  make masks and write
%   loop through all files ::
for iFile = filenames

    %%  make masks
    %   get filename ::
    filename = iFile{:};

    %   make mask ::
    [mask_data, variable, age, products] = make_cmip5_mask(filename, products, land_crit);

    %%  write data
    write_cmip5_variable(mask_data, variable, products, age, fullfile(exp_pro_path, age, strcat('land_mask_', age, '.nc')));
    
end

%%  end program
%%=========================================================================