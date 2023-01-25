% Take IS2 cross over locations & create a central location and define the strain diamond


clear variables 
close all
%% define everything needed to run code
plotting = 1; % if 1 = yes, if 0 = no
couplinglengthscale = 5; 
coordinatesign = -1; 
%% Read IS2 coordinates (variable centered_is2_locations = IS2 ROI = mean point of all laser crossovers; in

laser_xy = textread('is2locationdata/Pakistoq_Cycle_15_sample.txt');

cd ~/Documents/Projects/satellite_uplift/scripts/icestrainrates/

% Find the indicies of the points in a cluster
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

%pathtovelocity = 'icevelocitydata\'; % Christian's path
pathtovelocity = '~/Documents/Projects/satellite_uplift/scripts/icestrainrates/icevelocitydata/'; % Lauren's path
velocityfiletype = '*nc';

vel = readvelocitydata(velocitydatasource, pathtovelocity, velocityfiletype);

%clearvars -except vel centered_is2_locations

%% Read ice thickness data

thicknessdatasource = 1; %1 = bedmachine, 2 = ??

%pathtothickness = 'BedMachine\'; % Christian's path
pathtothickness = '~/Documents/Data/BedMachine/v5/'; % Lauren's path
thicknessfiletype = '*nc';

thick = readicethicknessdata(thicknessdatasource, pathtothickness, thicknessfiletype);

%% %% Extract a regional thicknesses ---> put this into function

xmax = max(centered_is2_locations(:,1), [], "all") + 20*1000; % do we want a smaller box still?
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

%% Extract regional velocities and indexes of crossovers ---> put this into function
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
vel.v_vel = vel.v_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
% For uncertainty, incorporate e_vel_std and n_vel_std here from readvelocitydata

tmpvelxy(:,1) = vel.x (:); 
tmpvelxy(:,2) = vel.y(:);

for ii = 1:length(centered_is2_locations)
    [index2(ii), distance2(ii)] = knnsearch(tmpvelxy(:,1:2), centered_is2_locations(ii,:), 'K', 1,'distance', 'euclidean');
    [vel_row(ii), vel_col(ii)] = ind2sub(size(vel.x),index2(ii)); 
end

vel_row = vel_row';
vel_col = vel_col';

clearvars -except vel centered_is2_locations thick thicknessesoflocations vel_row vel_col plotting couplinglengthscale laser_xy coordinatesign

%% Calculate strain rate

%crop the gridded data to each IS2 location

%% Setting Constants and Pre-allocating Matrices Outside The Loop

% Set time steps
dt=vel.tbands(2) - vel.tbands(1);
nt=dt * 2;%time step for advection and number of steps---CJT find the most appropriate nt this should also be in DAYS

% Set the rotation angle
alphas = 0:9:90; %rotates 9* every loop
eresults = -9999 + zeros(length(alphas), 3); % intialized using value not likely to occur

% Pre-allocate Matrices
% Temporary values to find nans
tmp_xp =   -9999*ones(length(centered_is2_locations),5,nt);
tmp_yp =   -9999*ones(length(centered_is2_locations),5,nt);
tmp_up =   -9999*ones(length(centered_is2_locations),5,nt); 
tmp_vp =   -9999*ones(length(centered_is2_locations),5,nt); 
% Saved values through rotation/advection loops
xpisave =  -9999*ones(length(alphas),5); 
ypisave =  -9999*ones(length(alphas),5);
lof1save = -9999*ones(length(alphas),1);
lof2save = -9999*ones(length(alphas),1);
xpsave =   -9999*ones(length(alphas),5);
ypsave =   -9999*ones(length(alphas),5);

for ii = 1 %:length(centered_is2_locations)

%% Set Constants Needed Within Loop
ROIhalflength = thicknessesoflocations(ii) .* couplinglengthscale;
gridcellcount = ceil(ROIhalflength ./ vel.pixelsize);  %decide if round, ceiling or floor is the best option

%% Renaming variables
% Define the X and Y grid of the velocity data
X = vel.x;
Y = vel.y;

% Define the velocity grids and fill velocity field NaNs
u = vel.e_vel;
v = vel.n_vel;
u_interp = inpaint_nans(u); % Can change method 1-8: https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
v_interp = inpaint_nans(v);

%% Define an inital bounding box for each IS2 ROI/crossover location before rotation using Vel_row and Ve_col index
bottom         = [vel_row(ii), vel_col(ii) + gridcellcount];
right       = [vel_row(ii) + gridcellcount, vel_col(ii)];
left        = [vel_row(ii) - gridcellcount, vel_col(ii)];
top      = [vel_row(ii), vel_col(ii) - gridcellcount];
centloc     = [vel_row(ii), vel_col(ii)];

%% Initial lengths based on initial coordinates
Xpi = [X(centloc(1), centloc(2)), X(top(1), top(2)),  X(right(1), right(2)), X(bottom(1), bottom(2)), X(left(1), left(2))];
Ypi = [Y(centloc(1), centloc(2)), Y(top(1), top(2)), Y(right(1), right(2)), Y(bottom(1), bottom(2)), Y(left(1), left(2))];

% Leg lengths of the intial strain diamond
a10=sqrt((Xpi(5)-Xpi(1))^2+(Ypi(5)-Ypi(1))^2);
a20=sqrt((Xpi(4)-Xpi(1))^2+(Ypi(4)-Ypi(1))^2);
b10=sqrt((Xpi(5)-Xpi(2))^2+(Ypi(5)-Ypi(2))^2);
b20=sqrt((Xpi(4)-Xpi(3))^2+(Ypi(4)-Ypi(3))^2);
c10=sqrt((Xpi(2)-Xpi(1))^2+(Ypi(2)-Ypi(1))^2);
c20=sqrt((Xpi(3)-Xpi(1))^2+(Ypi(3)-Ypi(1))^2);
d10=sqrt((Xpi(4)-Xpi(2))^2+(Ypi(4)-Ypi(2))^2);
d20=sqrt((Xpi(5)-Xpi(3))^2+(Ypi(5)-Ypi(3))^2);

halflength = a10; % l = a10 = a20 = c10 = c20;

for jj = 1:length(alphas)
    
    alpha = alphas(jj);
    
    lof1=halflength*cosd(alpha-90);
    lof2=halflength*sind(alpha-90);
    lof1save(jj,:) = lof1;
    lof2save(jj,:) = lof2;

    % X and Y now need to change with each rotation
    % Initial values before advection
    Xpi = [X(centloc(1), centloc(2)), (X(centloc(1), centloc(2)) -lof1), (X(centloc(1), centloc(2)) +lof1), (X(centloc(1), centloc(2)) +lof2), (X(centloc(1), centloc(2)) -lof2)];
    Ypi = [Y(centloc(1), centloc(2)), (Y(centloc(1), centloc(2)) -lof2), (Y(centloc(1), centloc(2)) +lof2), (Y(centloc(1), centloc(2)) -lof1), (Y(centloc(1), centloc(2)) +lof1)];
    Xp=Xpi;
    Yp=Ypi;
    xpisave(jj,:) = Xpi;
    ypisave(jj,:) = Ypi;
   
    for it=1:nt % Advecting for nt time-steps

        % Integrating dx/dt=u,dy/dt=v over a time step dt using Improved Euler Integration
        % up, vp - velocities at positions at beginning of time-step
        up=interp2(Y,X,u_interp,Yp,Xp);
        vp=interp2(Y,X,v_interp,Yp,Xp);

       !!%lauren deal with % Find NaN values in up/vp during advection
        tmp_up(ii,:,it) = interp2(Y,X,u,Yp,Xp);
        tmp_vp(ii,:,it) = interp2(Y,X,v,Yp,Xp);
        tmp_xp(ii,:,it) = Xp;
        tmp_yp(ii,:,it) = Yp;

        % First step of Improved Euler is Forward Euler (this is the final position of the velocity data for the entire time window)
        Xpstar=Xp+up*dt;
        Ypstar=Yp+vp*dt;

        % Second step of Improved Euler
        upstar=interp2(Y,X,u_interp,Ypstar,Xpstar);
        vpstar=interp2(Y,X,v_interp,Ypstar,Xpstar);

        % Final positions after each Improved Euler Time Step
        Xp=Xp+(up+upstar)/2*dt;
        Yp=Yp+(vp+vpstar)/2*dt;

    end % for it=1:nt%advecting for nt time-steps
   
    % Save Xp/Yp values for each rotation
    xpsave(jj,:) = Xp;
    ypsave(jj,:) = Yp;  

    % Final lengths based on final coordinates
    a1f=sqrt((Xp(5)-Xp(1))^2+(Yp(5)-Yp(1))^2);
    a2f=sqrt((Xp(4)-Xp(1))^2+(Yp(4)-Yp(1))^2);
    b1f=sqrt((Xp(5)-Xp(2))^2+(Yp(5)-Yp(2))^2);
    b2f=sqrt((Xp(4)-Xp(3))^2+(Yp(4)-Yp(3))^2);
    c1f=sqrt((Xp(2)-Xp(1))^2+(Yp(2)-Yp(1))^2);
    c2f=sqrt((Xp(3)-Xp(1))^2+(Yp(3)-Yp(1))^2);
    d1f=sqrt((Xp(4)-Xp(2))^2+(Yp(4)-Yp(2))^2);
    d2f=sqrt((Xp(5)-Xp(3))^2+(Yp(5)-Yp(3))^2);
    
    % Calculating logarithmic strain rates
    tott=nt*dt; % Total time for virtual advection
    edot0=0.5/tott*(log(a1f/a10)+log(a2f/a10));
    edot45=0.5/tott*(log(b1f/b10)+log(b2f/b10));
    edot90=0.5/tott*(log(c1f/c10)+log(c2f/c10));
    edot135=0.5/tott*(log(d1f/d10)+log(d2f/d10));

%% Final calculation of strain rate components
% a legs
ca=cosd(alpha);sa=sind(alpha);
% b legs
cam45=cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%cos(alpha-45)
sam45=-cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%sin(alpha-45)
% c legs
cam90=sind(alpha);%cos(alpha-90)
sam90=-cosd(alpha);%sin(alpha-90)
% d legs
cap45=cosd(alpha)/sqrt(2)-sind(alpha)/sqrt(2);%cos(alpha+45)
sap45=cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%sin(alpha+45)
%
A=[ca^2,2*ca*sa,sa^2;cam45^2,2*cam45*sam45,sam45^2;...
    cam90^2,2*cam90*sam90,sam90^2;cap45^2,2*cap45*sap45,sap45^2];
lsqmatrix=inv(A'*A)*A';
edotvector=lsqmatrix*[edot0;edot45;edot90;edot135];

% Components of the strain rate tensor
exxdot=edotvector(1);
exydot=edotvector(2);
eyydot=edotvector(3);
eresults(jj,:)=[exxdot,exydot,eyydot]; % Write estimated components into row of eresults

end %for jj = 1:length(alphas)

end %for ii = 1:length(centered_is2_locations)

%% Final Figures (If plotting multiple xover locations, move above outer for loop end)
% Figure 1
if plotting == 1
    figure
    hold on
    greenlandmap
    plot(xpisave(:,1,:),ypisave(:,1,:),'ko') %center
    plot(xpisave(:,2,:),ypisave(:,2,:),'mo') %top
    plot(xpisave(:,3,:),ypisave(:,3,:),'go') %right
    plot(xpisave(:,4,:),ypisave(:,4,:),'ro') %bottom
    plot(xpisave(:,5,:),ypisave(:,5,:),'co') %left
    plot(xpsave(:,1,:),ypsave(:,1,:),'k*') %center
    plot(xpsave(:,2,:),ypsave(:,2,:),'m*') %top
    plot(xpsave(:,3,:),ypsave(:,3,:),'g*') %right
    plot(xpsave(:,4,:),ypsave(:,4,:),'r*') %bottom
    plot(xpsave(:,5,:),ypsave(:,5,:),'c*') %left
    %
    plot(centered_is2_locations(:,1), centered_is2_locations(:,2),'b*')
    plot(Xpi(:), Ypi(:), '.r')
    plot(Xpstar(:), Ypstar(:), '.g')
    quiver(X(1:1:end, 1:1:end),Y(1:1:end, 1:1:end),u(1:1:end, 1:1:end),v(1:1:end, 1:1:end))
    streamline(X(1:1:end, 1:1:end),Y(1:1:end, 1:1:end),u(1:1:end, 1:1:end),v(1:1:end, 1:1:end))
    %
    % for tt = 1:length(xpisave)
    %     plot(xpisave(tt,:), ypisave(tt,:), '.', 'markersize', 15 )
    % end
    hold off
end %if plotting

%%
% Figure 2
if plotting == 1
    figure
    hold on
    %plot(Xpi, Ypi)
    a1 = scatter(xpisave(1,:), ypisave(1,:), 150, 1:5, "filled", 'square', 'markeredgecolor', 'b', 'LineWidth',3);
    scatter(xpisave(3,:), ypisave(3,:), 150, 1:5, "filled", 'markeredgecolor', 'b', 'LineWidth',3)
    scatter(xpisave(5,:), ypisave(5,:), 150, 1:5, "filled", 'markeredgecolor', 'b', 'LineWidth',3)
    scatter(xpisave(7,:), ypisave(7,:), 150, 1:5, "filled", 'markeredgecolor', 'b', 'LineWidth',3)
    scatter(xpisave(9,:), ypisave(9,:), 150, 1:5, "filled", 'markeredgecolor', 'b', 'LineWidth',3)
    scatter(xpisave(11,:), ypisave(11,:), 150, 1:5, "filled", 'diamond', 'markeredgecolor', 'b', 'LineWidth',3)

    b1 = scatter(xpsave(1,:), ypsave(1,:), 150, 1:5, "filled", 'square', 'markeredgecolor', 'r', 'LineWidth',3);
    scatter(xpsave(3,:), ypsave(3,:), 150, 1:5, "filled", 'markeredgecolor', 'r', 'LineWidth',3)
    scatter(xpsave(5,:), ypsave(5,:), 150, 1:5, "filled", 'markeredgecolor', 'r', 'LineWidth',3)
    scatter(xpsave(7,:), ypsave(7,:), 150, 1:5, "filled", 'markeredgecolor', 'r', 'LineWidth',3)
    scatter(xpsave(9,:), ypsave(9,:), 150, 1:5, "filled", 'markeredgecolor', 'r', 'LineWidth',3)
    scatter(xpsave(11,:), ypsave(11,:), 150, 1:5, "filled", 'diamond', 'markeredgecolor', 'r', 'LineWidth',3)

    !!%lauren, the strain diamonds seem more stretched using this set up, is this
    %due to using scatter?
    plot(centered_is2_locations(:,1), centered_is2_locations(:,2),'b*')
    quiver(X(1:1:end, 1:1:end),Y(1:1:end, 1:1:end),u(1:1:end, 1:1:end),v(1:1:end, 1:1:end))
    streamline(X(1:1:end, 1:1:end),Y(1:1:end, 1:1:end),u(1:1:end, 1:1:end),v(1:1:end, 1:1:end))

    legend([a1, b1],'Initial location', 'Advected/final location', 'location', 'southeast' )

    colormap(brewermap(5, 'spectral'))
    clim([1 6])
    cb2 = colorbar('Ticks',[1.5  2.5  3.5 4.5 5.5],...
        'TickLabels',{'Center','Top','Right','Bottom','Left'}, 'ticklength', 0);
end % if plotting

toc




%%
% xpi0 = [X(centloc(1), centloc(2)), (X(bottom(1), bottom(2)) ), (X(right(1), right(2)) ), (X(top(1), top(2)) ), (X(left(1), left(2)) )];
% ypi0 = [Y(centloc(1), centloc(2)), (Y(bottom(1), bottom(2)) ), (Y(right(1), right(2)) ), (Y(top(1), top(2)) ), (Y(left(1), left(2)))];
% 
% lof2save2 = (lof2save) + l;
% lof1save2 = (lof1save) + l;
% 
% for ii = 1:length(lof2save)
% 
% xpi1(ii,:) = [X(centloc(1), centloc(2)), (X(bottom(1), bottom(2))-lof1save(ii) ), (X(right(1), right(2))+lof1save(ii) ), (X(top(1), top(2))+lof2save(ii) ), (X(left(1), left(2))-lof2save(ii) )];
% ypi1(ii,:) = [Y(centloc(1), centloc(2)), (Y(bottom(1), bottom(2))-lof2save(ii) ), (Y(right(1), right(2))+lof2save(ii) ), (Y(top(1), top(2)) -lof1save(ii) ), (Y(left(1), left(2))+lof1save(ii))];
% 
% xpi2(ii,:) = [X(centloc(1), centloc(2)), (X(bottom(1), bottom(2))-lof1save2(ii) ), (X(right(1), right(2))+lof1save2(ii) ), (X(top(1), top(2))+lof2save2(ii) ), (X(left(1), left(2))-lof2save2(ii) )];
% ypi2(ii,:) = [Y(centloc(1), centloc(2)), (Y(bottom(1), bottom(2))-lof2save2(ii) ), (Y(right(1), right(2))+lof2save2(ii) ), (Y(top(1), top(2)) -lof1save2(ii) ), (Y(left(1), left(2))+lof1save2(ii))];
%   
% 
% end
% 
% figure 
% hold on 
% %plot(Xpi, Ypi)
% 
%  
% 
% scatter(xpi0(1,:), ypi0(1,:), 250, 1:5, "filled")
% scatter(xpi1(1,:), ypi1(1,:), 250, 1:5, "filled")
% scatter(xpi1(3,:), ypi1(3,:), 250, 1:5, "filled")
% scatter(xpi1(6,:), ypi1(6,:), 250, 1:5, "filled")
% 
% scatter(xpi2(1,:), ypi2(1,:), 250, 1:5, "filled", '*')
% scatter(xpi2(3,:), ypi2(3,:), 250, 1:5, "filled", '*')
% scatter(xpi2(6,:), ypi2(6,:), 250, 1:5, "filled", '*')
% 
% 
% colormap(brewermap(5, 'spectral'))
% caxis([1 6])
% cb2 = colorbar('Ticks',[1.5  2.5  3.5 4.5 5.5],...
%          'TickLabels',{'Center','Top','Right','Bottom','Left'}, 'ticklength', 0)
% 
% 
% %     figure
%     hold on
%     plot(Xpi, Ypi)
%     scatter(Xpi, Ypi, 250, 1:5, "filled")
%     colormap(brewermap(5, 'spectral'))
%     caxis([1 6])
%     cb2 = colorbar('Ticks',[1.5  2.5  3.5 4.5 5.5],...
%         'TickLabels',{'Center','Top','Right','Bottom','Left'}, 'ticklength', 0)
%     hold off
% % 

