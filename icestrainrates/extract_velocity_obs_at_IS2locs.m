% take IS2 cross over locations & create a central location and define the
% strain diamond
%lca 12-08-2022

clear variables 
close all
clc

%% Read IS2 coordinates

laser_xy = textread('is2locationdata/Pakistoq_Cycle_15_sample.txt');


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
        meanmatrix(ii,jj) = mean(laser_xy(tmp,jj),1, 'omitnan');
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

pathtovelocity = 'icevelocitydata';
velocityfiletype = '*nc';


v = readvelocitydata(velocitydatasource, pathtovelocity, velocityfiletype);

%% 
















