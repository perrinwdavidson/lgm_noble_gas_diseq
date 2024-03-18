function cmip_model_data = read_cmip5_variable(variable, products, paths)
%--------------------------------------------------------------------------
%   purpose: compiling multiple cmip data into a cell
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 06.07.22
%--------------------------------------------------------------------------
%%  let me know what is going o
disp(append('Aggregating raw CMIP5 ', upper(variable), ' model data from:')); 

%%  retrieve important stats
%   get number of models to read ::
NUMMOD = size(products, 1); 

%%  read data
%   pre-allocate final model array ::
cmip_model_data = cell(1, NUMMOD); 

%   loop through all models and read data ::
for iMod = 1 : 1 : NUMMOD

    %   write out how we are doing ::
	disp(products{iMod, 1});

	%    set path to directory ::
	iPath = fullfile(paths{1}, products{iMod, 1 : 2}, paths{2}); 

	%    get all files in directory ::
	filenames = dir(iPath);
	filenames = filenames(~[filenames.isdir], :);

	%   get number of files ::
	NUMFILES = length(filenames);

	%   run through all files and store in initial cell ::
    for iFile = 1 : 1 : NUMFILES

		%   make filename ::
		get_filename = fullfile(filenames(iFile).folder, filenames(iFile).name);

		%   get variable data ::
    	if iFile == 1

            %   get values ::
       		get_file.value = ncread(get_filename, variable);  % need to add in selection of pressure

            %   get shape, to decide if sic or ws ::
            file_shape = length(size(get_file.value));

            %   selection 100000 Pa (closest to sea level) if ws ::
            if file_shape == 4

                get_file.value = squeeze(get_file.value(:, :, 1, :));

            end

        else

            %   if sic, then just save ::
            if file_shape == 3

                get_file.value = cat(3, get_file.value, ncread(get_filename, variable));

            %   if ws, then need to grab sea level ::
            elseif file_shape == 4

                %   get the full array ::
                get_value = ncread(get_filename, variable);

                %   select only the sea level pressure ::
                get_file.value = cat(3, get_file.value, squeeze(get_value(:, :, 1, :)));

            end

        end

        %   show some files ::
        if mod(iFile, 10) == 0

            disp(append('Done reading file ', num2str(iFile)));

        end

    end

	%   get coordinates ::
	get_file.lon = double(ncread(get_filename, 'lon'));
   get_file.lat = double(ncread(get_filename, 'lat'));

	%   store data ::
	cmip_model_data{iMod} = get_file;

end

%--------------------------------------------------------------------------
end