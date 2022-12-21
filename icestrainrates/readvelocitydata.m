function v = readvelocitydata(datasource, path1, filetype)


if datasource == 1
    
    location1 = [path1, filetype];
    
    files        = dir(location1); %will load all the nc4 files in the folder
    files        = struct2cell(files);
    files        = files(1, 1:end);
    files        = files';

    files2 = [path1, files{1}];
    %This information is specific to promice data
    v.x      = ncread(files2, 'x');
    v.y      = ncread(files2, 'y');
    [v.x, v.y] = meshgrid(v.x, v.y); %these lines make 
    v.x = v.x'; v.y=v.y';

    for ii = 1:length(files)
        clear files2
        files2 = [path1, files{ii}];
        v.time(:,:,ii)   = ncread(files2, 'time') + datenum(1990,1,1);  %This may need to be changed for python
        v.tbands(:,:,ii) = ncread(files2, 'time_bnds') + datenum(1990,1,1); %This may need to be changed for python
        v.e_vel(:,:,ii)  = ncread(files2, 'land_ice_surface_easting_velocity');
        v.n_vel(:,:,ii)  = ncread(files2, 'land_ice_surface_northing_velocity');
        v.v_vel(:,:,ii)  = ncread(files2, 'land_ice_surface_vertical_velocity');

    end


elseif datasource == 2
%

elseif datasource == 3



else 
    disp('The code does not understand your choice of velocity data: Abort!')

end 