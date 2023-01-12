% take IS2 cross over locations & create a central location and define the
% strain diamond
%lca 12-08-2022, edited lca & CJT 01/11/23

clear variables 
close all
%% define everything needed to run code
plotting = 1; % if 1 = yes, if 0 = no
couplinglengthscale = 1; 
%% Read IS2 coordinates (variables centered_is3_locations = IS2 ROI = mean point of all laser crossovers; in

laser_xy = textread('is2locationdata/Pakistoq_Cycle_15_sample.txt');

cd ~/Documents/Projects/satellite_uplift/scripts/icestrainrates/

% find the indicies of the points in a cluster
for ii = 1:length(laser_xy)
%[nearestneighbors, D] = knnsearch(distance(ii,ii), distance(ii,:)', 'K', 4,'distance', 'euclidean')
[index, distance] = knnsearch(laser_xy/1000, laser_xy(ii,:)/1000, 'K', 4,'distance', 'euclidean');
  
for jj = 1:4
    if distance(jj) > 1
        
       distance(jj) = 0;
       index(jj) = 0;
    end
   end

indexmatrix(ii,:) = index;
%distancematrix(ii,:)    = distance;

end

% we might have to sort the distances too! if the distances are needed

sortedindexmatrix = sort(indexmatrix,2);

sortedindexmatrix= unique(sortedindexmatrix, 'rows');
sortedindexmatrix(sortedindexmatrix==0) = NaN;


for ii = 1:length(sortedindexmatrix)
    tmp = sortedindexmatrix(ii,:);
    tmp = tmp(~isnan(tmp));
    for jj = 1:2 %if the height is added, this becomes 3
        centered_is2_locations(ii,jj) = mean(laser_xy(tmp,jj),1, 'omitnan');
    end
end

%clearvars -except centered_is2_locations   



%% Read velocity data (variable vel = northing and easting vel data and gridded x and y polar stereo) 

velocitydatasource = 1; % 1 = promice velocity ; 2 = go live velocity (give the actual data name) ; 3 a fun new dataset we dont know about yet

% pathtovelocity = 'icevelocitydata\'; % Christian's path
pathtovelocity = '~/Documents/Projects/satellite_uplift/scripts/icestrainrates/icevelocitydata/'; % Lauren's path
velocityfiletype = '*nc';

vel = readvelocitydata(velocitydatasource, pathtovelocity, velocityfiletype);

%clearvars -except vel centered_is2_locations

%% read ice thickness data

thicknessdatasource = 1; %1 = bedmachine, 2 = ??

%pathtothickness = 'BedMachine\'; % Christian's path
pathtothickness = '~/Documents/Data/BedMachine/v5/'; % Lauren's path
thicknessfiletype = '*nc';

thick = readicethicknessdata(thicknessdatasource, pathtothickness, thicknessfiletype);


%% %% extract a regional thicknesses ---> put this into function

xmax = max(centered_is2_locations(:,1), [], "all") + 20*1000; 
xmin = min(centered_is2_locations(:,1), [], "all") - 20*1000;
ymax = max(centered_is2_locations(:,2), [], "all") + 20*1000; 
ymin = min(centered_is2_locations(:,2), [], "all") - 20*1000;

tmpthickx = double(thick.x);
tmpthicky = double(thick.y);
  tmpthickx((tmpthickx<xmin)) = NaN;
  tmpthickx((tmpthickx>xmax)) = NaN;
  tmpthicky((tmpthicky<ymin)) = NaN;
  tmpthicky((tmpthicky>ymax)) = NaN;

tmponest = tmpthickx + tmpthicky;
tmponest(~isnan(tmponest)) = 1;

tmpthick.x = tmpthickx .* tmponest;
tmpthick.x(any(isnan(tmpthick.x), 2), :) = [];
tmpthick.y = tmpthicky .* tmponest;
tmpthick.y(any(isnan(tmpthick.y), 2), :) = [];
tmpthick.t = thick.thickness .* tmponest;
tmpthick.t(any(isnan(tmpthick.t), 2), :) = [];

tmpthickxy(:,1) = tmpthickx(:); 
tmpthickxy(:,2) = tmpthicky(:);
tmpthickxy(:,3) = thick.thickness(:);
tmpthickxy(any(isnan(tmpthickxy), 2), :) = [];

for ii = 1:length(centered_is2_locations)
    [index2(ii), distance2(ii)] = knnsearch(tmpthickxy(:,1:2), centered_is2_locations(ii,:), 'K', 1,'distance', 'euclidean');
end

thicknessesoflocations = tmpthickxy(index2,3); %in METERS


%% extract a regional velocities and indexes of crossovers ---> put this into function
tmpvelx = double(vel.x);
tmpvely = double(vel.y);
  tmpvelx((tmpvelx<xmin)) = NaN;
  tmpvelx((tmpvelx>xmax)) = NaN;
  tmpvely((tmpvely<ymin)) = NaN;
  tmpvely((tmpvely>ymax)) = NaN;

tmponev = tmpvelx + tmpvely;
tmponev(~isnan(tmponev)) = 1;

% CJT --> include the other velocity variables here
vel.x = vel.x .* tmponev;
indicies1 = find(tmponev(:)==1);
[~, maxxind] = max(indicies1, [], 'omitnan');
[~, minxind] = min(indicies1, [], 'omitnan');
 [upleft(1,1), upleft(1,2)] = ind2sub(size(tmponev),minxind); 
 [botright(1,1), botright(1,2)] = ind2sub(size(tmponev),maxxind); 


%vel.x(any(isnan(vel.x), 2), :) = [];
vel.y = vel.y .* tmponev;
%vel.y(any(isnan(vel.y), 2), :) = [];
vel.n_vel = vel.n_vel .* tmponev;
vel.n_vel(any(isnan(vel.n_vel), 2), :) = []; %ERROR - Matrix index is out of range for deletion.
vel.e_vel = vel.e_vel .* tmponev;
vel.e_vel(any(isnan(vel.e_vel), 2), :) = [];

tmpvelxy(:,1) = tmpvelx(:); 
tmpvelxy(:,2) = tmpvely(:);
% tmpvelxy(:,3) = t.thickness(:);
% tmpvelxy(any(isnan(tmpthickxy), 2), :) = [];

for ii = 1:length(centered_is2_locations)
    [index2(ii), distance2(ii)] = knnsearch(tmpvelxy(:,1:2), centered_is2_locations(ii,:), 'K', 1,'distance', 'euclidean');
    [vel_row(ii), vel_col(ii)] = ind2sub(size(tmponev),index2(ii)); 
end

vel_row = vel_row';
vel_col = vel_col';
%clearvars -except vel centered_is2_locations thick thicknessesoflocations vel_row vel_col plotting couplinglengthscale laser_xy

%% Use hari's code to calculate strain rate

%crop the gridded data to each IS2 location

for ii = 1 %:length(centered_is2_locations)

%define a bounding box for each IS2 ROI/crossover location

ROIhalflength = thicknessesoflocations(ii) .* couplinglengthscale;

gridcellcount = ceil(ROIhalflength ./ vel.pixelsize);  %decide if round, ceiling or floor is the best option

indlocCent  = [vel_row(ii), vel_col(ii)];
%these labels might not be correct depending on the sign of x and y (e.g.
%the sign of the northing and easting in the projected velocity data, but
%the grid will still center correctly
topleft     = [vel_row(ii) - gridcellcount, vel_col(ii) + gridcellcount]
topright    = [vel_row(ii) + gridcellcount, vel_col(ii) + gridcellcount]
bottomleft  = [vel_row(ii) - gridcellcount, vel_col(ii) - gridcellcount]
bottomright = [vel_row(ii) + gridcellcount, vel_col(ii) - gridcellcount]

%%
if plotting ==1
    figure
    hold on
    %greenlandmap
    surface(vel.x, vel.y, vel.v_vel)
    scatter(laser_xy(:,1), laser_xy(:,2))
%     plot( centered_is2_locations(1,1), centered_is2_locations(1,2)  , 'g*')
%     plot(X(:),Y(:), '.g')
%     plot(vel.x(topleft(1), topleft(2)),vel.y(topleft(1), topleft(2)), 'r*')
%     plot(vel.x(topright(1), topright(2)),vel.y(topright(1), topright(2)), 'b*')
%     plot(vel.x(bottomleft(1), bottomleft(2)),vel.y(bottomleft(1), bottomleft(2)), 'g*')
%     plot(vel.x(bottomright(1), bottomright(2)),vel.y(bottomright(1), bottomright(2)), 'm*')
%     plot(vel.x(indlocCent(1), indlocCent(2)),vel.y(indlocCent(1), indlocCent(2)), 'pm')
%     %axis([-1.37e5 -1.30e5 -2.27e6 -2.255e6])
end
%%

X = vel.x(topleft(1):bottomright(1), bottomright(2):topleft(2));
Y = vel.y(topleft(1):bottomright(1), bottomright(2): topleft(2));










% define the X and Y grid of the velocity data
%------making a hypotherical grid and hypothetical velocity dataset
%------hypothetical x-y grid: %x=0:0.01:1;y=0:0.01:1;
X = vel.x;
Y = vel.y;

% define the velocity grids
u = vel.e_vel;
v = vel.n_vel;

%% Figure to check the field
if plotting == 1
    % CJT edit as needed.
    figure
    hold on
    greenlandmap
    quiver(X(1:10:end, 1:10:end),Y(1:10:end, 1:10:end),u(1:10:end, 1:10:end),v(1:10:end, 1:10:end))
    streamline(X(1:10:end, 1:10:end),Y(1:10:end, 1:10:end),u(1:10:end, 1:10:end),v(1:10:end, 1:10:end))

end 

%% 




end %for ii = 1:length(centered_is2_locations)








