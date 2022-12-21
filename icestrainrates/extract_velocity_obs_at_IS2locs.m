% take IS2 cross over locations & create a central location and define the
% strain diamond
%lca 12-08-2022

clear variables 
close all
%clc

%% Read IS2 coordinates

laser_xy = textread('is2locationdata/Pakistoq_Cycle_15_sample.txt');

cd ~/Documents/Projects/satellite_uplift/scripts/icestrainrates/

%% Quick check with a plot
% 
% figure 
% hold on
% greenlandmap
% scatter(points(:,1), points(:,2))

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
tmpthickx(:,1) = double(t.x);
tmpthicky(:,1) = double(t.y);
  tmpthickx((tmpthickx<xmin)) = NaN;
  tmpthickx((tmpthickx>xmax)) = NaN;
  tmpthicky((tmpthicky<ymin)) = NaN;
  tmpthicky((tmpthicky>ymax)) = NaN;


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


tmpvelxy(:,1) = tmpthickx(:);
tmpvelxy(:,2) = tmpthicky(:);
for ii = 1:length(v)
    tmpthickxy(:,ii+2) = v.e_vel;
end



for ii = 1:length(v)




end


%% figure to check

figure 
hold on
greenlandmap
scatter(tmpthickxy(:,1), tmpthickxy(:,2), 40, tmpthickxy(:,3), 'filled', 'MarkerEdgeColor', 'none')

scatter(laser_xy(:,1), laser_xy(:,2))
scatter(centered_is2_locations(:,1), centered_is2_locations(:,2), 40, 'r*')
scatter([xmin, xmax, xmax, xmin], [ymax, ymax, ymin, ymin], 60, 'g^')














