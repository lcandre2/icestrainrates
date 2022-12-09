%read PRMOICE velocity data and calculate the strain rates 

clear variables 
close all
clc

%velocity is in meters per day
%%

cd ~/Documents/Data/IceVel_Promice/test/

files        = dir('*.nc'); %will load all the nc4 files in the folder
files        = struct2cell(files);
files        = files(1, 1:end);
files   = files';

v.x      = ncread(files{1}, 'x');
v.y      = ncread(files{1}, 'y');
[v.x, v.y] = meshgrid(v.x, v.y);
v.x = v.x'; v.y=v.y';

for ii = 1:length(files)

    v.time(:,:,ii)   = ncread(files{ii}, 'time') + datenum(1990,1,1);
    v.tbands(:,:,ii) = ncread(files{ii}, 'time_bnds') + datenum(1990,1,1);
    v.e_vel(:,:,ii)  = ncread(files{ii}, 'land_ice_surface_easting_velocity');
    v.n_vel(:,:,ii)  = ncread(files{ii}, 'land_ice_surface_northing_velocity');
    v.v_vel(:,:,ii)  = ncread(files{ii}, 'land_ice_surface_vertical_velocity');

end




%%



%% calculate strain rates
cd /Users/lcandre2/Documents/Projects/satellite_uplift/scripts

for ii = 1:length(files)

 disp(['Starting strain rate calculations for... ', datestr(v.time(:,:,ii))] )

 [v.elon(:,:,ii), v.etrans(:,:,ii), v.eshear(:,:,ii), v.eEff(:,:,ii), v.ez(:,:,ii)] ...
     = calculate_strainrates(v.e_vel(:,:,ii), v.n_vel(:,:,ii), NaN, 500, 10^-4, 1, ...
                             500, 500, 8, NaN, NaN,NaN);
 disp(['Done with... ', datestr(v.time(:,:,ii))] )

end

%%


for kk = 1:length(files)
    figure('position', [ 1          536        1306         441])
    subplot(1,3,1)
    hold on
    greenlandmap
    surf(v.x, v.y, 365.25*sqrt(v.e_vel(:,:,kk).^2 +  v.n_vel(:,:,kk).^2)        , 'edgecolor', 'none')
    cb = colorbar
    cmocean('haline', 30)
    caxis([0 300])
    ylabel(cb, 'Ice speed m y^{-1}')


    subplot(1,3,2)
    hold on
    greenlandmap
    surf(v.x, v.y, v.elon(:,:,kk)*365.25         , 'edgecolor', 'none')
    title([datestr(v.tbands(1,1,kk)), ' to ', datestr(v.tbands(2,1,kk)) ]  )
    cb = colorbar
    caxis([-0.08 0.08])
    cmocean('balance', 20, 'pivot',0)
   
     ylabel(cb, 'Longitudinal strain rate y^{-1}')

    subplot(1,3,3)
    hold on
    greenlandmap
    surf(v.x, v.y, v.ez(:,:,kk)*365.25         , 'edgecolor', 'none')
    cb = colorbar
    caxis([-0.08 0.08])
    cmocean('curl', 20, 'pivot',0)
    ylabel(cb, 'Vertical strain rate y^{-1}')
end



%% plot histograms

for mm = 1:length(files)
    clear veltmp lontmp vtmp
    veltmp = 365.25*sqrt(v.e_vel(:,:,mm).^2 +  v.n_vel(:,:,mm).^2);
    veltmp = veltmp(:);

    lontmp = v.elon(:,:,mm)*365.25;
    lontmp = lontmp(:);

    vtmp = v.ez(:,:,mm)*365.25;
    vtmp = vtmp(:);

    figure('position', [ 1         462        1179         515])
    subplot(1,3,1)
    hold on
    histogram(log(veltmp), 'edgecolor', 'none')
    ylabel('Count')
    xlabel('Log(ice velocity)')

    subplot(1,3,2)
    hold on
    histogram(lontmp, 'edgecolor', 'none')
    axis([-0.2 0.2 0 10^6])
    ylabel('Count')
    xlabel('Longitudinal strain rate')
    title([datestr(v.tbands(1,1,mm)), ' to ', datestr(v.tbands(2,1,mm)) ]  )


    subplot(1,3,3)
    hold on
    histogram(vtmp, 'edgecolor', 'none')
    axis([-0.2 0.2 0 10^6])
    ylabel('Count')
    xlabel('Vertical strain rate')
end





