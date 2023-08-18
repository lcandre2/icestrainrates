function [v, veloc_sum] = readvelocitydata(datasource, path1, filetype, is2roi)


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
        v.tot_vel(:,:,ii)  = ncread(files2, 'land_ice_surface_velocity_magnitude');
        v.vert_vel(:,:,ii)  = ncread(files2, 'land_ice_surface_vertical_velocity');
        v.tot_vel_std(:,:,ii)  = ncread(files2, 'land_ice_surface_velocity_magnitude_std');
        v.e_vel_std(:,:,ii)  = ncread(files2, 'land_ice_surface_easting_velocity_std');
        v.n_vel_std(:,:,ii)  = ncread(files2, 'land_ice_surface_northing_velocity_std');
        v.e_vel_base(:,:,ii)  = ncread(files2, 'land_ice_surface_easting_velocity');
        v.n_vel_base(:,:,ii)  = ncread(files2, 'land_ice_surface_northing_velocity');

    end


v.pixelsize = abs(v.x(1) - v.x(2)); %This makes the implicit assumption that the velocity data have 'square' pixels! 


elseif datasource == 2
%
        location1 = [path1, filetype];

        files        = dir(location1); %will load all the nc4 files in the folder
        files        = struct2cell(files);
        files        = files(1, 1:end);
        files        = files';

        %IS2 Crossover dTs
        % xover_dT_file = 'is2locationdata\Crossover_dTs\dT_Petermann_Cycle_2_crossovers.txt';
         
        
        for ii = 1:length(files)
            % if files(ii)(18:25)
            clear files2
            files2 = [path1, files{ii}];
            
            x      = ncread(files2, 'x');
            y      = ncread(files2, 'y');
            
            for jj = 1:length(is2roi)
                if and(is2roi(jj,1) >min(x), is2roi(jj,1)<max(x)) && and(is2roi(jj,2) >min(y), is2roi(jj,2)<max(y) )
                    disp('inside');
                    veloc(ii) = 1;
                else
                    disp('outside')
                    veloc(ii)  = 0;
                end

            end
            
            veloc_sum = sum(veloc);

            [v{ii}.x, v{ii}.y] = meshgrid(x, y); %these lines make
            v{ii}.x = v{ii}.x'; v{ii}.y =v{ii}.y';
            v{ii}.x(:,:,ii) = v{ii}.x; v{ii}.y(:,:,ii) = v{ii}.y;

            date_center = ncreadatt(files2,'img_pair_info','date_center');
            v{ii}.date_center(:,:,ii) = date_center(1:8);
            v{ii}.date_dt(:,:,ii) = ncreadatt(files2,'img_pair_info','date_dt');
            v{ii}.start_time(:,:,ii) = datetime(v{ii}.date_center(1:8),'Format','uuuuMMdd') - (v{ii}.date_dt)/2;
            v{ii}.end_time(:,:,ii)   = datetime(v{ii}.date_center(1:8), 'Format','uuuuMMdd') + (v{ii}.date_dt)/2;
           
            
            v{ii}.tot_vel(:,:,ii) = ncread(files2,'v');
            v{ii}.e_vel(:,:,ii) = ncread(files2,'vx');
            v{ii}.n_vel(:,:,ii) = ncread(files2,'vy');
            v{ii}.tot_vel_std(:,:,ii) = ncread(files2,'v_error');
            e_vel_std = ncreadatt(files2,'vx','error');
            n_vel_std = ncreadatt(files2,'vy','error');
            v{ii}.e_vel_std(:,:,ii) = ones(height(v{ii}.x),width(v{ii}.y)) * e_vel_std;
            v{ii}.n_vel_std(:,:,ii) = ones(height(v{ii}.x),width(v{ii}.y)) * n_vel_std;
            v{ii}.e_vel_base(:,:,ii) = ncread(files2,'vx');
            v{ii}.n_vel_base(:,:,ii) = ncread(files2,'vy');

        end

        %Get average of all files in ii?


v{ii}.pixelsize(:,:,ii) = abs(v{ii}.x(1) - v{ii}.x(2)); %This makes the implicit assumption that the velocity data have 'square' pixels! 


elseif datasource == 3



else 
    disp('The code does not understand your choice of velocity data: Abort!')

end 