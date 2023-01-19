% take IS2 cross over locations & create a central location and define the
% strain diamond
%lca 12-08-2022, edited lca & CJT 01/11/23, edited CJT 01/16/23

clear variables 
close all
%% define everything needed to run code
plotting = 1; % if 1 = yes, if 0 = no
couplinglengthscale = 5; 
%% Read IS2 coordinates (variables centered_is3_locations = IS2 ROI = mean point of all laser crossovers; in

laser_xy = textread('is2locationdata/Pakistoq_Cycle_15_sample.txt');

%cd ~/Documents/Projects/satellite_uplift/scripts/icestrainrates/

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

pathtovelocity = 'icevelocitydata\'; % Christian's path
%pathtovelocity = '~/Documents/Projects/satellite_uplift/scripts/icestrainrates/icevelocitydata/'; % Lauren's path
velocityfiletype = '*nc';

vel = readvelocitydata(velocitydatasource, pathtovelocity, velocityfiletype);

%clearvars -except vel centered_is2_locations

%% read ice thickness data

thicknessdatasource = 1; %1 = bedmachine, 2 = ??

pathtothickness = 'BedMachine\'; % Christian's path
%pathtothickness = '~/Documents/Data/BedMachine/v5/'; % Lauren's path
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

indicies1 = find(tmponev(:)==1);
[upleft(1,1), upleft(1,2)] = ind2sub(size(tmponev),indicies1(1)) ;
[botright(1,1), botright(1,2)] = ind2sub(size(tmponev),indicies1(end));

vel.y = vel.y(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.x = vel.x(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.n_vel = vel.n_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.e_vel = vel.e_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.vert_vel = vel.vert_vel(upleft(1):botright(1), upleft(2): botright(2)) ; %
%CJT update readvelocitydata for 132


tmpvelxy(:,1) = vel.x (:); 
tmpvelxy(:,2) = vel.y(:);

for ii = 1:length(centered_is2_locations)
    [index2(ii), distance2(ii)] = knnsearch(tmpvelxy(:,1:2), centered_is2_locations(ii,:), 'K', 1,'distance', 'euclidean');
    [vel_row(ii), vel_col(ii)] = ind2sub(size(vel.x),index2(ii)); 
end

vel_row = vel_row';
vel_col = vel_col';



clearvars -except vel centered_is2_locations thick thicknessesoflocations vel_row vel_col plotting couplinglengthscale laser_xy

%% Use hari's code to calculate strain rate

%crop the gridded data to each IS2 location

for ii = 1 %:length(centered_is2_locations)

%%
% define the X and Y grid of the velocity data
%------making a hypotherical grid and hypothetical velocity dataset
%------hypothetical x-y grid: %x=0:0.01:1;y=0:0.01:1;
X = vel.x;
Y = vel.y;

% define the velocity grids
u = vel.e_vel;
v = vel.n_vel;

%%


%define a bounding box for each IS2 ROI/crossover location

ROIhalflength = thicknessesoflocations(ii) .* couplinglengthscale;

gridcellcount = ceil(ROIhalflength ./ vel.pixelsize);  %decide if round, ceiling or floor is the best option


%these labels might not be correct depending on the sign of x and y (e.g.
%the sign of the northing and easting in the projected velocity data, but
%the grid will still center correctly
topleft     = [vel_row(ii) - gridcellcount, vel_col(ii) + gridcellcount];
topright    = [vel_row(ii) + gridcellcount, vel_col(ii) + gridcellcount];
bottomleft  = [vel_row(ii) - gridcellcount, vel_col(ii) - gridcellcount];
bottomright = [vel_row(ii) + gridcellcount, vel_col(ii) - gridcellcount];
centloc  = [vel_row(ii), vel_col(ii)];

% fix to set up diamond to be not the corners, but the top middle, left
% middle, etc.
%%CJT check to see if lines 180-184 are correctly defining the strain
%%diamond


%% CJT check to see if in the right order on lines 188 and 189
Xpi = [X(centloc(1), centloc(2)), X(topleft(1), topleft(2)), X(topright(1), topright(2)), X(bottomright(1), bottomright(2)), X(bottomleft(1), bottomleft(2))];
Ypi = [Y(centloc(1), centloc(2)), Y(topleft(1), topleft(2)), Y(topright(1), topright(2)), Y(bottomright(1), bottomright(2)), Y(bottomleft(1), bottomleft(2))];


Xp=Xpi;
Yp=Ypi;

% if plotting == 1
figure
hold on
plot(Xp,Yp,'r*')

%____________
%CJT dt and nt need to be defined based on the velocity data
%_____________

%currently dt is in DAYS, just like the velocity data is in m/DAY.
dt=vel.tbands(2) - vel.tbands(1);

% Set to 6 hours, *2 = 12
nt=dt * 2;%time step for advection and number of steps---CJT find the most appropriate nt this should also be in DAYS

% Uses inpaint_nans to interpolate and extrapolate NaN elements in 2d array (u and v)
%u_interp = inpaint_nans(u); % Can change method 1-8: https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
%v_interp = inpaint_nans(v);

% Creates a 4 column vector if X, Y, u, and v for scatteredInterpolant 
q2(:,1) = X(:);q2(:,2) = Y(:);q2(:,3) = u(:);q2(:,4) = v(:);

for it=1:nt%advecting for nt time-steps

    %integrating dx/dt=u,dy/dt=v over a time step dt using Improved Euler Integration
    Fu1 = scatteredInterpolant(q2(:,1),q2(:,2),q2(:,3),"linear");
    Fv1 = scatteredInterpolant(q2(:,1),q2(:,2),q2(:,4),"linear");
    up = Fu1(Xp,Yp);
    vp = Fv1(Xp,Yp);
   
    % First step of Improved Euler is Forward Euler (this is the final position of the velocity data for the entire time window)
    Xpstar=Xp+up*dt;
    Ypstar=Yp+vp*dt;

    % Second step of Improved Euler
    Fu = scatteredInterpolant(q2(:,1),q2(:,2),q2(:,3),"linear");
    Fv = scatteredInterpolant(q2(:,1),q2(:,2),q2(:,4),"linear");
    upstar = Fu(Xpstar,Ypstar);
    vpstar = Fv(Xpstar,Ypstar);

    %final positions after each Improved Euler Time Step
    Xp=Xp+(up+upstar)/2*dt;
    Yp=Yp+(vp+vpstar)/2*dt;

    % ERROR PROBLEM: ScatteredInterpolant is supposed to fill in NaN
    % values but it is not doing so. NaNs remain - meaning the use of
    % inpaint_nans would be required regardless of this vs. interp2. (Lines
    % 217 and 218
end
%

%plot final position after nt timesteps
% 
% if plotting == 1
%     % CJT edit as needed.
%     figure
%     hold on
    greenlandmap
    plot(Xpi(:), Ypi(:), '.r')
    plot(Xpstar(:), Ypstar(:), '.g')
%     
    quiver(X(1:1:end, 1:1:end),Y(1:1:end, 1:1:end),u(1:1:end, 1:1:end),v(1:1:end, 1:1:end))
    streamline(X(1:1:end, 1:1:end),Y(1:1:end, 1:1:end),u(1:1:end, 1:1:end),v(1:1:end, 1:1:end))
    plot(Xp,Yp,'ro')
% end 
% 
% 
%now on to the strain rate calculations - calculating the lengths of
%different segments in the strain diamond (a,b,c,d)
%initial lengths based on initial coordinates
a10=sqrt((Xpi(5)-Xpi(1))^2+(Ypi(5)-Ypi(1))^2);
a20=sqrt((Xpi(4)-Xpi(1))^2+(Ypi(4)-Ypi(1))^2);
b10=sqrt((Xpi(5)-Xpi(2))^2+(Ypi(5)-Ypi(2))^2);
b20=sqrt((Xpi(4)-Xpi(3))^2+(Ypi(4)-Ypi(3))^2);
c10=sqrt((Xpi(2)-Xpi(1))^2+(Ypi(2)-Ypi(1))^2);
c20=sqrt((Xpi(3)-Xpi(1))^2+(Ypi(3)-Ypi(1))^2);
d10=sqrt((Xpi(4)-Xpi(2))^2+(Ypi(4)-Ypi(2))^2);
d20=sqrt((Xpi(5)-Xpi(3))^2+(Ypi(5)-Ypi(3))^2);
%final lengths based on final coordinates
a1f=sqrt((Xp(5)-Xp(1))^2+(Yp(5)-Yp(1))^2);
a2f=sqrt((Xp(4)-Xp(1))^2+(Yp(4)-Yp(1))^2);
b1f=sqrt((Xp(5)-Xp(2))^2+(Yp(5)-Yp(2))^2);
b2f=sqrt((Xp(4)-Xp(3))^2+(Yp(4)-Yp(3))^2);
c1f=sqrt((Xp(2)-Xp(1))^2+(Yp(2)-Yp(1))^2);
c2f=sqrt((Xp(3)-Xp(1))^2+(Yp(3)-Yp(1))^2);
d1f=sqrt((Xp(4)-Xp(2))^2+(Yp(4)-Yp(2))^2);
d2f=sqrt((Xp(5)-Xp(3))^2+(Yp(5)-Yp(3))^2);
%next calculating logarithmic strain rates
tott=nt*dt;%total time for virtual advection
edot0=0.5/tott*(log(a1f/a10)+log(a2f/a10));
edot45=0.5/tott*(log(b1f/b10)+log(b2f/b10));
edot90=0.5/tott*(log(c1f/c10)+log(c2f/c10));
edot135=0.5/tott*(log(d1f/d10)+log(d2f/d10));
%finally calculate strain rate components
lsqmatrix=[-1/4,1/4,3/4,1/4;0,1/2,0,1/2;3/4,1/4,-1/4,1/4];
edotvector=lsqmatrix*[edot0;edot45;edot90;edot135];
%components of the strain rate tensor
exxdot=edotvector(1)
exydot=edotvector(2)
eyydot=edotvector(3)
hold off


end %for ii = 1:length(centered_is2_locations)






