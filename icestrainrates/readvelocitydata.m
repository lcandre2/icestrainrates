function v = readvelocitydata(datasource, path, filetype)


if datasource == 1
    
    cd(path)
    
    files        = dir(filetype); %will load all the nc4 files in the folder
    files        = struct2cell(files);
    files        = files(1, 1:end);
    files        = files';


    %This information is specific to promice data
    v.x      = ncread(files{1}, 'x');
    v.y      = ncread(files{1}, 'y');
    [v.x, v.y] = meshgrid(v.x, v.y);
    v.x = v.x'; v.y=v.y';

    for ii = 1:length(files)

        v.time(:,:,ii)   = ncread(files{ii}, 'time') + datenum(1990,1,1);  %This may need to be changed for python
        v.tbands(:,:,ii) = ncread(files{ii}, 'time_bnds') + datenum(1990,1,1); %This may need to be changed for python
        v.e_vel(:,:,ii)  = ncread(files{ii}, 'land_ice_surface_easting_velocity');
        v.n_vel(:,:,ii)  = ncread(files{ii}, 'land_ice_surface_northing_velocity');
        v.v_vel(:,:,ii)  = ncread(files{ii}, 'land_ice_surface_vertical_velocity');

    end


elseif velocitydatasource == 2
%

elseif velocitydatasource == 3



else 
    disp('The code does not understand your choice of velocity data: Abort!')

end 