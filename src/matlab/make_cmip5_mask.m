function [mask_data, variable, age, products] = make_cmip5_mask(filename, products, land_crit)
%--------------------------------------------------------------------------
%   purpose: making cmip5 and uvic masks.
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 22.07.22
%--------------------------------------------------------------------------
%%  configure
%   get variable name ::
variable = filename(7 : 11);

%   get age ::
age = filename(13 : 15); 

%   let me know what is going on ::
disp(append('Making land masks for ', upper(age)));

%%  load uvic data
%   load ::
load(fullfile('data', 'exp_raw', 'uvic', 'grid'), 'x', 'y', 'ideep'); 

%   make grid ::
[uvic_lat, uvic_lon] = meshgrid(y, x); 

%%  get inputs for this model
%   get variables names ::
group_names = ncread(filename, 'group_names');
variable_names = ncread(filename, 'variable_names');

%   append uvic ::
group_names = [group_names; {append('UVic ', upper(age), ' Default')}];
variable_names = [variable_names; {append('UVic Default ', upper(age), ' ' , upper(variable))}];

%   modify products to include uvic ::
products = [products; {'UVic', 'Default'}];

%   get number of models ::
NUMMOD = size(group_names, 1);

%   set number of months ::
NUMMON = 12;
        
%%  preallocate arrays
mask_data = cell([1, NUMMOD]); 

%%  loop through all models
%   make masks ::
for iMod = 1 : 1 : NUMMOD

    %   for cmip5 ::
    if iMod < 8

        %   get model data ::
        model_data = ncread(filename, append(group_names{iMod}, variable_names{iMod}));
        model_lon = ncread(filename, append(group_names{iMod}, 'lon'));
        model_lat = ncread(filename, append(group_names{iMod}, 'lat'));

        %   normalize to [0 1] ::
        if sum(model_data > 1, 'all') > 0

            model_data = model_data ./ 100;  % assume out of 100

        end

        %   make mask ::
        model_data_mask = model_data; 
        model_data(model_data_mask <= land_crit) = 1; 
        model_data(model_data_mask > land_crit) = 0; 

    %   for uvic ::
    elseif iMod == 8

        %   make mask ::
        ideep(ideep ~= 0) = 1; 
        ind_lm = find(ideep(:) == 1);
        land_mask = logical(ideep(:));
        clear('ideep');

        %   get model data ::
        model_data = land_mask;
        model_lon = uvic_lon; 
        model_lat = uvic_lat;

    end

    %   save data ::
    mask_data{iMod}.value = model_data;
    mask_data{iMod}.lon = model_lon; 
    mask_data{iMod}.lat = model_lat; 

    %   write out how we are doing ::
    disp(append('Done making mask for ', variable_names{iMod, 1}));

end

%--------------------------------------------------------------------------
end