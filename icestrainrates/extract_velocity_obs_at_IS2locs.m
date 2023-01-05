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

%% extract a regional thicknesses ---> put this into function

xmax = max(centered_is2_locations(:,1), [], "all") + 20*1000; 
xmin = min(centered_is2_locations(:,1), [], "all") - 20*1000;
ymax = max(centered_is2_locations(:,2), [], "all") + 20*1000; 
ymin = min(centered_is2_locations(:,2), [], "all") - 20*1000;

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

thicknessesoflocations = tmpthickxy(index2,3); %in METERS


%% extract a regional velocities and indexes of crossovers ---> put this into function
tmpvelx = double(v.x);
tmpvely = double(v.y);
  tmpvelx((tmpvelx<xmin)) = NaN;
  tmpvelx((tmpvelx>xmax)) = NaN;
  tmpvely((tmpvely<ymin)) = NaN;
  tmpvely((tmpvely>ymax)) = NaN;

tmponev = tmpvelx + tmpvely;
tmponev(~isnan(tmponev)) = 1;

tmpvel.x = tmpvelx .* tmponev;
tmpvel.x(any(isnan(tmpvel.x), 2), :) = [];
tmpvel.y = tmpvely .* tmponev;
tmpvel.y(any(isnan(tmpvel.y), 2), :) = [];
tmpvel.vy = v.n_vel .* tmponev;
tmpvel.vy(any(isnan(tmpvel.vy), 2), :) = [];
tmpvel.vx = v.e_vel .* tmponev;
tmpvel.vx(any(isnan(tmpvel.vx), 2), :) = [];

tmpvelxy(:,1) = tmpthickx(:); 
tmpvelxy(:,2) = tmpthicky(:);
% tmpvelxy(:,3) = t.thickness(:);
% tmpvelxy(any(isnan(tmpthickxy), 2), :) = [];

for ii = 1:length(centered_is2_locations)
    [index2(ii), distance2(ii)] = knnsearch(tmpvelxy(:,1:2), centered_is2_locations(ii,:), 'K', 1,'distance', 'euclidean');
    [row(ii), col(ii)] = ind2sub(size(tmponev),index2(ii)); %little bit more work here...
end


%thicknessesoflocations = tmpthickxy(index2,3); %in METERS




%% Now extract immediate thickness data 
% provide a loose bounding box to increase computational efficiency (this
% just increases the bounding box to 10km around each cross over location

couplinglengthscale = 1; % The number of ice thicknesses over which to calculate the r. FOR a 4 loop, this becomes 0.5:0.5:10
pixel_size = v.pixelsize; %in METERS
locMult = 2;

for ii = 1 %this for loop is to eventually do multiple cross over regions
    maxR = thicknessesoflocations * 10000 ./ pixel_size; %the max expected coupling lengthscale??
    %boxsize = 10*1000; %km * 100m to give meters... bounding box, could also be set by ice thickness
    %boxsize = thicknessesoflocations(ii);

%     locxmax = centered_is2_locations(ii,1) + boxsize;
%     locxmin = centered_is2_locations(ii,1) - boxsize;
%     locymax = centered_is2_locations(ii,2) + boxsize;
%     locymin = centered_is2_locations(ii,2) - boxsize;
     
    length_scale = thicknessesoflocations * couplinglengthscale;
        r = round(length_scale/pixel_size); % Finds the nearest number of pixels to the given length scale;
        if r == 0 % If the length scale rounds to 0, set it to 1
            r = 1;
        elseif r > maxR
            r = maxR;
        end
        %r = cast(r,'double');
        %rGrid(i,j) = r;
        
        % Define the lengths of the segments at the beginning of each calculation
        l0a1 = r;
        l0a2 = r;
        l0b1 = r*sqrt(2);
        l0b2 = r*sqrt(2);
        l0c1 = r;
        l0c2 = r;
        l0d1 = r*sqrt(2);
        l0d2 = r*sqrt(2);

        % Local square dimensions
        locDim = (2*locMult*r)+1;
        locCent = ceil(locDim/2);

        % Assign local coordinates to the stakes around each strain square
        rowCoords = [locCent,locCent-r,locCent,locCent+r,locCent];
        colCoords = [locCent,locCent,locCent+r,locCent,locCent-r];
    
 %% Initialize calculations for each center pixel
        % Set the initial strain experienced by each strain segment to zero
        stota1 = 0;
        stota2 = 0;
        stotb1 = 0;
        stotb2 = 0;
        stotc1 = 0;
        stotc2 = 0;
        stotd1 = 0;
        stotd2 = 0;
        
        % Set the current length of each strain segment to the original
        % lengths comprising the strain square
        lLasta1 = l0a1;
        lLasta2 = l0a2;
        lLastb1 = l0b1;
        lLastb2 = l0b2;
        lLastc1 = l0c1;
        lLastc2 = l0c2;
        lLastd1 = l0d1;
        lLastd2 = l0d2;
        
        % Set the current rows and columns to the coordinates of the strain
        % square
        curRows = rowCoords;
        curCols = colCoords;

         % Extract an array around the center point (i,j) that represents
        % twice the dimensions of the strain square in order to make later
        % calculations. This is the "local grid"
        %i = row index of crossover point on velocity grid
        %j = column index of crossover point on velocity grid
        sqVx =v.e_vel((row(ii)-(locMult*r)):(row(ii)+(locMult*r)),(col(ii)-(locMult*r)):(col(ii)+(locMult*r)));
        sqVy = v.n_vel((row(ii)-(locMult*r)):(row(ii)+(locMult*r)),(col(ii)-(locMult*r)):(col(ii)+(locMult*r)));
        sqThick = thick((row(ii)-(locMult*r)):(row(ii)+(locMult*r)),(col(ii)-(locMult*r)):(col(ii)+(locMult*r)));

        % Extract an array around the center point (i,j) that represents
        % just the strain square in order to calculate the average velocity
        % at the center point and determine a reasonable time interval
        sqVxmean = vx((row(ii)-r):(row(ii)+r),(col(ii)-r):(col(ii)+r));
        sqVymean = vy((row(ii)-r):(row(ii)+r),(col(ii)-r):(col(ii)+r));
        [sqRows,sqCols] = size(sqVx);
% Calculate mean velocity
        meanX = nanmean(nanmean(sqVxmean));
        meanY = nanmean(nanmean(sqVymean));
        meanVel = sqrt(meanX^2+meanY^2);
        
        % Let the stakes move by approximately one tenth of the length scale
        time = 0.1*r*pixel_size/meanVel;
        time = min(time, time_max);

        dtOrig = pixel_size/meanVel*.05; % Initialize time step as the time it takes to move 
        % a twentieth of the pixel length, according to the average velocity
        dt = dtOrig;
        t = 0; % Initialize time tracker
        
        % Extract the x- and y-velocities at the original stake positions
        for k = 1:5
            curXVels(k) = sqVx(rowCoords(k),colCoords(k));
            curYVels(k) = sqVy(rowCoords(k),colCoords(k));
        end
        if any(isnan([curXVels,curYVels]))==1
            continue
            % Go to the next pixel in the for-loop if any of the current
            % velocities are NaNs
        end
        
        %% Let the stakes move and measure strain rates
            while t<time
                % Move the stakes according to the x- and y-velocities and the
                % time step. Calculate the new column and row positions.
                % Also calculate the column and row positions according to
                % the improved Euler method to check for accuracy.
                t = t+dt;
                for k = 1:5
                    newRowCoords(k) = curRows(k) - ydir*curYVels(k)*dt/pixel_size;
                    newColCoords(k) = curCols(k) + curXVels(k)*dt/pixel_size;
                    [checkXVels(k),checkYVels(k)] = locInterp2(curRows(k),curCols(k),sqVx,sqVy);
                    checkRowCoords(k) = curRows(k) - ydir*0.5*(curYVels(k)+checkYVels(k))*dt/pixel_size;
                    checkColCoords(k) = curCols(k) + 0.5*(curXVels(k)+checkXVels(k))*dt/pixel_size;
                    errorCriteriaY(k) = abs((checkRowCoords(k)-newRowCoords(k))/checkRowCoords(k));
                    errorCriteriaX(k) = abs((checkColCoords(k)-newColCoords(k))/checkColCoords(k));
                end
                if any(isnan([checkXVels,checkYVels]))==1
                    break
                end
                %% Adaptive time stepping
                if any([errorCriteriaY,errorCriteriaX] >= tol)==1
                    t = t-dt; % Reverse the time to what it was before
                    dt = dt/2; % Make a smaller time step
                    % We leave the current rows, columns, and velocities the
                    % same    
                else
                %% Make final calculations
                % Check to be sure that the stakes haven't moved outside of the
                % local grid
                    if max(newRowCoords) <= sqRows && max(newColCoords) <= sqCols && min(newRowCoords)>=1 && min(newColCoords)>=1
                        % Calculate the current length 
                        lfa1 = sqrt((newRowCoords(1) - newRowCoords(2))^2+(newColCoords(1)-newColCoords(2))^2);
                        lfa2 = sqrt((newRowCoords(1) - newRowCoords(4))^2+(newColCoords(1)-newColCoords(4))^2);
                        lfb1 = sqrt((newRowCoords(2) - newRowCoords(5))^2+(newColCoords(2)-newColCoords(5))^2);
                        lfb2 = sqrt((newRowCoords(3) - newRowCoords(4))^2+(newColCoords(3)-newColCoords(4))^2);
                        lfc1 = sqrt((newRowCoords(1) - newRowCoords(5))^2+(newColCoords(1)-newColCoords(5))^2);
                        lfc2 = sqrt((newRowCoords(1) - newRowCoords(3))^2+(newColCoords(1)-newColCoords(3))^2);
                        lfd1 = sqrt((newRowCoords(2) - newRowCoords(3))^2+(newColCoords(2)-newColCoords(3))^2);
                        lfd2 = sqrt((newRowCoords(4) - newRowCoords(5))^2+(newColCoords(4)-newColCoords(5))^2);

                        % Calculate the current strains and strain rates
                        stra1 = log(lfa1/lLasta1);
                        stra2 = log(lfa2/lLasta2);
                        strb1 = log(lfb1/lLastb1);
                        strb2 = log(lfb2/lLastb2);
                        strc1 = log(lfc1/lLastc1);
                        strc2 = log(lfc2/lLastc2);
                        strd1 = log(lfd1/lLastd1);
                        strd2 = log(lfd2/lLastd2);

                        ea1 = stra1/dt;
                        ea2 = stra2/dt;
                        eb1 = strb1/dt;
                        eb2 = strb2/dt;
                        ec1 = strc1/dt;
                        ec2 = strc2/dt;
                        ed1 = strd1/dt;
                        ed2 = strd2/dt;
                  
                        % Update the new rows and columns as current
                        curRows = newRowCoords;
                        curCols = newColCoords;
                        
                        % Update the current lengths as the previous
                        % lengths
                        lLasta1 = lfa1;
                        lLasta2 = lfa2;
                        lLastb1 = lfb1;
                        lLastb2 = lfb2;
                        lLastc1 = lfc1;
                        lLastc2 = lfc2;
                        lLastd1 = lfd1;
                        lLastd2 = lfd2;
                        
                        % Update the running total of strain for each
                        % segment
                        stota1 = stota1 + stra1;
                        stota2 = stota2 + stra2;
                        stotb1 = stotb1 + strb1;
                        stotb2 = stotb2 + strb2;
                        stotc1 = stotc1 + strc1;
                        stotc2 = stotc2 + strc2;
                        stotd1 = stotd1 + strd1;
                        stotd2 = stotd2 + strd2;
                        
                        % Calculate final strain rate components
                        ea1f = stota1/t;
                        ea2f = stota2/t;
                        eb1f = stotb1/t;
                        eb2f = stotb2/t;
                        ec1f = stotc1/t;
                        ec2f = stotc2/t;
                        ed1f = stotd1/t;
                        ed2f = stotd2/t;

                        % Reset the time step
                        dt = dtOrig;

                        % Set the current velocities to those calculated at
                        % the end of the time step
                        curXVels = checkXVels;
                        curYVels = checkYVels;

             
                    else
                        break
                    end
                    % If a stake has moved outside the strain square, leave the
                    % current rows and columns and current velocities the same
                    % as the previous time step, and simply move on to
                    % calculating the strain rate components. This will be the
                    % final value for this strain square. Change the time to
                    % kick it out of the while loop
                end
            end
            if t<.5*time % Don't calculate values if the stakes have been 
                % allowed to move for less than half of the designated time
                % increment
                exGrid(i,j) = nan;
                eyGrid(i,j) = nan;
                exyGrid(i,j) = nan;
                centerAlphas(i,j) = nan;
            else
            %% Create strain rate grids

            % Average the strain rate components
            ea = (ea1f + ea2f)/2;
            eb = (eb1f + eb2f)/2;
            ec = (ec1f + ec2f)/2;
            ed = (ed1f + ed2f)/2;
            
            % Calculate coordinate-oriented strain values
            exGrid(row(ii),col(ii)) = .25*(eb + ed - ea) + .75*ec;
            exyGrid(row(ii),col(ii)) = .5*eb - .5*ed;
            eyGrid(row(ii),col(ii)) = .75*ea + .25*(eb + ed - ec);
    

%% Calculate flow orientation
            % Calculate a grid of flow directions so that the grid-oriented
            % strain rates can be rotated outside of the for-loop to align with
            % local flow directions
            centerVelX = vx(row(ii),col(ii));
            centerVelY = vy(row(ii),col(ii));

            if centerVelX>0 && centerVelY>0
                centerAlphas(i,j) = atand(centerVelY/centerVelX);
            elseif centerVelX<0 && centerVelY>0
                centerAlphas(i,j) = atand(centerVelY/centerVelX)+180;
            elseif centerVelX<0 && centerVelY<0
                centerAlphas(i,j) = atand(centerVelY/centerVelX)+180;
            elseif centerVelX>0 && centerVelY<0
                centerAlphas(i,j) = atand(centerVelY/centerVelX)+360;
            elseif centerVelX>0 && centerVelY==0
                centerAlphas(i,j) = 0;
            elseif centerVelX==0 && centerVelY>0
                centerAlphas(i,j) = 90;
            elseif centerVelX<0 && centerVelY==0
                centerAlphas(i,j) = 180;
            elseif centerVelX==0 && centerVelY<0
                centerAlphas(i,j) = 270;
            end
            end
    end  



% %% figure to check
% 
% figure 
% hold on
% greenlandmap
% scatter(laser_xy(:,1), laser_xy(:,2))
% scatter(centered_is2_locations(:,1), centered_is2_locations(:,2), 40, 'r*')
% scatter([xmin, xmax, xmax, xmin], [ymax, ymax, ymin, ymin], 60, 'g^')
% 

% 
% 
% %% extract regional velocities
% 
% % only need to crop the first one with values
% tmpvelx = double(v.x);
% tmpvely = double(v.y);
% tmpvelx((tmpvelx<xmin)) = NaN;
% tmpvelx((tmpvelx>xmax)) = NaN;
% tmpvely((tmpvely<ymin)) = NaN;
% tmpvely((tmpvely>ymax)) = NaN;
% 
% % make a fake matrix of ones and nans
% 
% tmpones  = tmpvelx + tmpvely;
% tmpones(~isnan(tmpones)) = 1;
% 
% tmpvel.x = tmpvelx .* tmpones;
% tmpvel.y = tmpvely .* tmpones;
% tmpvel.x(any(isnan(tmpvel.x), 2), :) = [];
% tmpvel.y(any(isnan(tmpvel.y), 2), :) = [];
% 
% for ii = 1:length(v)
%     clear v1 v2 v3
%     v1 = v.e_vel .* tmpones;
%     %v1(any(isnan(v1), 2), :) = [];
% 
%     v2 = v.n_vel .* tmpones;
%     %v2(any(isnan(v2), 2), :) = [];
% 
%     v3 = v.v_vel .* tmpones;
%     %v3(any(isnan(v3), 2), :) = [];
% 
%     tmpvel.evel(:,:,ii) = v1;
%     tmpvel.nvel(:,:,ii) = v2;
%     tmpvel.vel(:,:, ii) = v3;
% end
% 
% %% 
% 
% figure 
% hold on
% greenlandmap
% scatter(tmpvelx(:), tmpvely(:), 10, tmpvel.vel(:), 'filled', 'markeredgecolor', 'none')
% scatter(laser_xy(:,1), laser_xy(:,2))
% scatter(centered_is2_locations(:,1), centered_is2_locations(:,2), 40, 'r*')
% scatter([xmin, xmax, xmax, xmin], [ymax, ymax, ymin, ymin], 60, 'g^')
% 
% 
% 
% %% Determine local flow lines 
% %1 must extend velocity field to coastline.
% 
% clearvars -except xmin xmax ymin ymax
% 
% 
% 
% %%
% % %% Calculate the longitudinal strain rate with varying length scales
% % 
% % edot = 1/l0 .* deltaL /deltaT;
% % 
% % %1. velocity data is in nothing easting - we need information in the
% % %direction of ice flow. so we need to reorient the n and e velocities to
% % %provide longitudinal velocity. This calculation will have a certain
% % %uncertainity (UC1).
% % %2. Calculate the longitudinal, lateral and vertical  strain rate over ice
% % %thickness. Lengthscale matters, generally use 2-10x ice thickness (UC2).
% % %3. the height change due to vertical strain is depedent on ice thickness (UC3). 
% % % uncertainty in the measured ice thickness, but the bigger uncertainty 
% % % is what ice thickness from the defined lengthscale to use. There is also
% % % uncertainity in the assumption the Ezz is uniform with depth.
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 









%%

% for ii = 1:length(v)
% tic
% 
% 
%     [v.elon(:,:,ii), v.etrans(:,:,ii), v.eshear(:,:,ii), v.eEff(:,:,ii), v.ez(:,:,ii)] ...
%      = calculatestrainrates(v.e_vel(:,:,ii), v.n_vel(:,:,ii), NaN, 500, 10^-4, 1, ...
%                              500, 500, 8, NaN, NaN,NaN);
% 
%     toc
% end
% 
% %% figure to check
% % 
% 
% 
% figure 
% hold on
% greenlandmap
%  scatter(tmpvel.x(:), tmpvel.y(:), 40, s.elon(:), 'filled', 'MarkerEdgeColor', 'none')
% 
% scatter(laser_xy(:,1), laser_xy(:,2))
% scatter(centered_is2_locations(:,1), centered_is2_locations(:,2), 40, 'r*')
% scatter([xmin, xmax, xmax, xmin], [ymax, ymax, ymin, ymin], 60, 'g^')
% 
% 
% 











