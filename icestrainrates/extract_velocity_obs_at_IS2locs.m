% take IS2 cross over locations & create a central location and define the
% strain diamond
%lca 12-08-2022

clear variables 
close all
%clc

%% Read IS2 coordinates

laser_xy = textread('is2locationdata/Pakistoq_Cycle_15_sample.txt');

cd ~/Documents/Projects/satellite_uplift/scripts/icestrainrates/



%% What we need to do.



%1. come up with mean/central value for all lasers
%2. velocities around the points 
    %two types of velocity data: full greenland or tiled (option 1: just
    %read all the files ; option 2: use coordinates from IS2 to subselect
    %the tiles of interest.
%3. ice thicknesses around points
%4. some way to define a length scale

%% Determine mean IS2 points 
% for ii = 1 :length(laser_xy)
% %for jj = 1:length(laser_xy) %the length here may be the max number of points in a cluster
%  distance(ii,:)  =sqrt( ((laser_xy(ii,1) - laser_xy(:,1)).^2) +  ((laser_xy(ii,2) - laser_xy(:,2)).^2 ) ) /1000;
% end 
%% find the indicies of the points in a cluster
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



%%
% figure 
% hold on
% %greenlandmap
% scatter(laser_xy(:,1), laser_xy(:,2))
% scatter(meanmatrix(:,1), meanmatrix(:,2), 40, 'r*')

%% Read velocity data

velocitydatasource = 1; % 1 = promice velocity ; 2 = go live velocity (give the actual data name) ; 3 a fun new dataset we dont know about yet

pathtovelocity = '~/Documents/Projects/satellite_uplift/scripts/icestrainrates/icevelocitydata/';
velocityfiletype = '*nc';


v = readvelocitydata(velocitydatasource, pathtovelocity, velocityfiletype);

%% read ice thickness data

thicknessdatasource = 1; %1 = bedmachine, 2 = ??

pathtothickness = '~/Documents/Data/BedMachine/v5/';
thicknessfiletype = '*nc';

t = readicethicknessdata(thicknessdatasource, pathtothickness, thicknessfiletype);

%% Now extract immediate thickness data 

% provide a loose bounding box to increase computational efficiency (this
% just increases the bounding box to 20km greater than the min and max of
% the is2 cross overs
xmax = max(centered_is2_locations(:,1), [], "all") + 20*1000; 
xmin = min(centered_is2_locations(:,1), [], "all") - 20*1000;
ymax = max(centered_is2_locations(:,2), [], "all") + 20*1000; 
ymin = min(centered_is2_locations(:,2), [], "all") - 20*1000;

%% figure to check

% figure 
% hold on
% greenlandmap
% scatter(laser_xy(:,1), laser_xy(:,2))
% scatter(centered_is2_locations(:,1), centered_is2_locations(:,2), 40, 'r*')
% scatter([xmin, xmax, xmax, xmin], [ymax, ymax, ymin, ymin], 60, 'g^')

%% extract an regional thicknesses
tmpthickx = double(t.x);
tmpthicky = double(t.y);
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
tmpthick.t = t.thickness .* tmponest;
tmpthick.t(any(isnan(tmpthick.t), 2), :) = [];

tmpthickxy(:,1) = tmpthickx(:); 
tmpthickxy(:,2) = tmpthicky(:);
tmpthickxy(:,3) = t.thickness(:);
tmpthickxy(any(isnan(tmpthickxy), 2), :) = [];

for ii = 1:length(centered_is2_locations)
    [index2(ii), distance2(ii)] = knnsearch(tmpthickxy(:,1:2), centered_is2_locations(ii,:), 'K', 1,'distance', 'euclidean');
end

thicknessesoflocations = tmpthickxy(index2,3);


%% extract regional velocities

% only need to crop the first one with values
tmpvelx = double(v.x);
tmpvely = double(v.y);
tmpvelx((tmpvelx<xmin)) = NaN;
tmpvelx((tmpvelx>xmax)) = NaN;
tmpvely((tmpvely<ymin)) = NaN;
tmpvely((tmpvely>ymax)) = NaN;

% make a fake matrix of ones and nans

tmpones  = tmpvelx + tmpvely;
tmpones(~isnan(tmpones)) = 1;

tmpvel.x = tmpvelx .* tmpones;
tmpvel.y = tmpvely .* tmpones;
tmpthick.x(any(isnan(tmpthick.x), 2), :) = [];
tmpthick.y(any(isnan(tmpthick.y), 2), :) = [];

for ii = 1:length(v)
    clear v1 v2 v3
    v1 = v.e_vel .* tmpones;
    v1(any(isnan(v1), 2), :) = [];

    v2 = v.n_vel .* tmpones;
    v2(any(isnan(v2), 2), :) = [];

    v3 = v.v_vel .* tmpones;
    v3(any(isnan(v3), 2), :) = [];

    tmpvel.evel(:,:,ii) = v1;
    tmpvel.nvel(:,:,ii) = v2;
    tmpvel.vel(:,:, ii) = v3;
end

%% Calculate the longitudinal strain rate with varying length scales

edot = 1/l0 .* deltaL /deltaT;

%1. velocity data is in nothing easting - we need information in the
%direction of ice flow. so we need to reorient the n and e velocities to
%provide longitudinal velocity. This calculation will have a certain
%uncertainity (UC1).
%2. Calculate the longitudinal, lateral and vertical  strain rate over ice
%thickness. Lengthscale matters, generally use 2-10x ice thickness (UC2).
%3. the height change due to vertical strain is depedent on ice thickness (UC3). 
% uncertainty in the measured ice thickness, but the bigger uncertainty 
% is what ice thickness from the defined lengthscale to use. There is also
% uncertainity in the assumption the Ezz is uniform with depth.






















%%

for ii = 1:length(v)
tic


    [v.elon(:,:,ii), v.etrans(:,:,ii), v.eshear(:,:,ii), v.eEff(:,:,ii), v.ez(:,:,ii)] ...
     = calculatestrainrates(v.e_vel(:,:,ii), v.n_vel(:,:,ii), NaN, 500, 10^-4, 1, ...
                             500, 500, 8, NaN, NaN,NaN);

    toc
end

%% figure to check
% 


figure 
hold on
greenlandmap
 scatter(tmpvel.x(:), tmpvel.y(:), 40, s.elon(:), 'filled', 'MarkerEdgeColor', 'none')

scatter(laser_xy(:,1), laser_xy(:,2))
scatter(centered_is2_locations(:,1), centered_is2_locations(:,2), 40, 'r*')
scatter([xmin, xmax, xmax, xmin], [ymax, ymax, ymin, ymin], 60, 'g^')














