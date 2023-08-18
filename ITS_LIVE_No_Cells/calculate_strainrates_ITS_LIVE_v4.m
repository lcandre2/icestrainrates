% calculate strain rates for IceSat-2 

clear variables
%close all

%'Should' only need to modify through line 52!!

%% define inputs
%user:
islauren = 0; %use lauren's paths if set to 1 , otherwise, use Christian's paths

%% script features (CHECK CHECK CHECK)

couplinglengthscale = (1:3);%(1:0.5:5); %determines how many and what lengthscales strain rates are calculated for.

velocityerrrorquant = 'montecarlo'; %'none'; 'std'; 'montecarlo' %this tells the script what type of velocity error analysis to use
numsim = 8; %100 % Velocity Perturbations - This is only used for the monte carlo simulations (the results will always be numsim +1 since the first calc will always be for the reported velocity
numstd = 1; %use either 1 standard deviation or two. This applies to 'std' and 'montecarlo' error analysis. This needs to be aligned with the reported velocity errors.
% 

largeboundingboxsize = 30; %km, crops the vel and thick data to be more reasonable

saveoutputs = 1; % if 1 = yes, if 0 = no -- Save output .mat files to the current directory for each centered_is2_location

%% IS2 inputs

Cycle = 2; % Change to correct IS2 cycle number based on laser_xy input for filename purposes
Region = 'Petermann_Terminus'; % Change region or glacier name for filename purposes


if islauren == 1
    is2pathname = '~/Documents/Projects/satellite_uplift/scripts/icestrainrates/is2locationdata/';
else
    is2pathname = 'is2locationdata\'; %Pakistoq_Cycle_15_sample.txt';
end
is2filetype = '*.txt';

%% velocity inputs

velocitydatasource = 2; % 1 = promice velocity ; 2 = ITS_LIVE velocity ; 3 a fun new dataset we dont know about yet

if velocitydatasource == 1
    velocityfiletype = '*nc'; % Promice
elseif velocitydatasource == 2
    velocityfiletype = '*nc'; % ITS_LIVE
end

if islauren == 1
    pathtovelocity = '~/Documents/Projects/satellite_uplift/scripts/icestrainrates/icevelocitydata/'; % Lauren's path
else
    pathtovelocity = 'icevelocitydata\';%Petermann-ITS_LIVE-c02\';
end
   
%% thickness inputs

thicknessdatasource = 1; %1 = bedmachine, 2 = ??
thicknessfiletype = '*nc';

if islauren == 1
    pathtothickness = '~/Documents/Data/BedMachine/v5/'; % Lauren's path
else
    pathtothickness =  'BedMachine\'; 
end

%% Calculate the centered is2 cross over location (ROI) and load relevant velocity and thickness data
%These are all functions to read the data

is2roi = calculateis2crossoverroi(is2pathname, is2filetype); %This function is not well tested. Check if files are different than using during testing.

[vel,veloc_sum] = readvelocitydata(velocitydatasource, pathtovelocity, velocityfiletype, is2roi);
    if veloc_sum == 0
        disp('Stopping script due to mismatched IS2 and velocity coordinates...');
        return
    end

thick  = readicethicknessdata(thicknessdatasource, pathtothickness, thicknessfiletype);

% vel.e_vel_base = vel.e_vel;
% vel.n_vel_base = vel.n_vel;

%% crop velocity data and extract thickness and vel coordinates at IS2 ROI

%This is quite slow...
[roithick, roibederror, velrow, velcol, vel] = is2roiinfo(is2roi, thick, vel, largeboundingboxsize);


%% assign velocity uncertainties

    if strcmp(velocityerrrorquant, 'none')

        disp('No error analysis chosen, code will proceed')

    elseif strcmp(velocityerrrorquant, 'std')

        vel.e_vel(:,:,1) = vel.e_vel + vel.e_vel_std.*numstd;
        vel.e_vel(:,:,2) = vel.e_vel - vel.e_vel_std.*numstd;
        vel.n_vel(:,:,1) = vel.n_vel + vel.n_vel_std.*numstd;
        vel.n_vel(:,:,2) = vel.n_vel - vel.n_vel_std.*numstd;

    elseif strcmp(velocityerrrorquant, 'montecarlo')

        %this elseif statement creates a bivariate normal distribution of
        %velocities using the number of simulations (numsim) and

        size1 = [size(vel.x), numsim];
        test  = randn(size1);
        test2 = randn(size1);

        tmpe(:,:,1) = vel.e_vel;
        tmpe(:,:,2) = vel.e_vel + vel.e_vel_std.*numstd;
        tmpe(:,:,3) = vel.e_vel - vel.e_vel_std.*numstd;
        tmpn(:,:,1) = vel.n_vel;
        tmpn(:,:,2) = vel.n_vel + vel.n_vel_std.*numstd;
        tmpn(:,:,3) = vel.n_vel - vel.n_vel_std.*numstd;

        edist = test.*vel.e_vel_std;
        ndist = test2.*vel.n_vel_std; % independent of each other, not true but conservative
        %edist = (tmpe(:,:,2) - tmpe(:,:,3)).*test + tmpe(:,:,3); %scaled randn values to fall within range of +- 1 std
        %ndist = (tmpn(:,:,2) - tmpn(:,:,3)).*test2 + tmpn(:,:,3);

        vel.e_vel = vel.e_vel_base + edist;
        vel.n_vel = vel.n_vel_base + ndist;
        %need to keep parent velocity field
        vel.e_vel = cat(3,vel.e_vel, edist);
        vel.n_vel = cat(3, vel.n_vel, ndist);

        clear tmpe tmpn edist ndist test test2 size1


    else
        disp('You need to determine what velocity error structure you would like to use!!!!!')
        return
    end
disp('Finished Assigning Velocity Error Values (Line 136), Moving on to Strain Rate Calcuations')


%% Calculate strain rates
% Set time steps
% strainrates = struct('vp',[],'crossoverROI',[],'crossoverthick',[],'roibederror',[]);

if velocitydatasource == 1
    dt=vel.tbands(2) - vel.tbands(1);
    nt=dt * 4; % Time step for advection and number of steps---CJT find the most appropriate nt this should also be in DAYS

elseif velocitydatasource == 2
    dt=vel.date_dt;
    nt=dt/2;% * 4; % Time step for advection and number of steps in days.
end
 
% Set the rotation angle
alphas = 0:9:90; %rotates 9* every loop
%eresults = -9999 + zeros(length(alphas), 3); % intialized using value not likely to occur

X = vel.x;
Y = vel.y;

for ii = 1:2%length(is2roi)%:13 %length(is2roi) %This is for each crossover location
    tic
    for jj = 1:numsim+1 %this is for each velocity perturbation

        u = vel.e_vel(:,:,jj);
        v = vel.n_vel(:,:,jj); %v = vel.e_vel(:,:,jj);

        u_interp = inpaint_nans(u); % Can change method 1-8: https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
        v_interp = inpaint_nans(v);

        for  kk = 1:length(couplinglengthscale)

            %define the coupling lengthscale based on the grid
            gridcellcount =...
                            ceil((roithick(ii) .* couplinglengthscale(kk) )...
                            ./ vel.pixelsize);


            if gridcellcount > largeboundingboxsize
                disp('Coupling length is larger than the bounding box, Abort')
                return
            end

            % Define an inital bounding box for each IS2 ROI/crossover location before rotation using Vel_row and Ve_col index
            bottom   = [velrow(ii), velcol(ii) + gridcellcount];
            right    = [velrow(ii) + gridcellcount, velcol(ii)];
            left     = [velrow(ii) - gridcellcount, velcol(ii)];
            top      = [velrow(ii), velcol(ii) - gridcellcount];
            centloc  = [velrow(ii), velcol(ii)];

            Xpi      = [X(centloc(1), centloc(2)), X(top(1), top(2)),...
                X(right(1), right(2)), X(bottom(1), bottom(2)), X(left(1), left(2))];
            Ypi      = [Y(centloc(1), centloc(2)), Y(top(1), top(2)),...
                Y(right(1), right(2)), Y(bottom(1), bottom(2)), Y(left(1), left(2))];

            a10      = sqrt((Xpi(5)-Xpi(1))^2+(Ypi(5)-Ypi(1))^2);
            a20      = sqrt((Xpi(4)-Xpi(1))^2+(Ypi(4)-Ypi(1))^2);
            b10      = sqrt((Xpi(5)-Xpi(2))^2+(Ypi(5)-Ypi(2))^2);
            b20      = sqrt((Xpi(4)-Xpi(3))^2+(Ypi(4)-Ypi(3))^2);
            c10      = sqrt((Xpi(2)-Xpi(1))^2+(Ypi(2)-Ypi(1))^2);
            c20      = sqrt((Xpi(3)-Xpi(1))^2+(Ypi(3)-Ypi(1))^2);
            d10      = sqrt((Xpi(4)-Xpi(2))^2+(Ypi(4)-Ypi(2))^2);
            d20      = sqrt((Xpi(5)-Xpi(3))^2+(Ypi(5)-Ypi(3))^2);

            halflength = a10; % l = a10 = a20 = c10 = c20;

            for mm = 1:length(alphas)

                alpha   = alphas(mm);

                lof1    = halflength*cosd(alpha-90);
                lof2    = halflength*sind(alpha-90);

                Xpi     = [X(centloc(1), centloc(2)), (X(centloc(1),...
                    centloc(2)) -lof1), (X(centloc(1), centloc(2)) +lof1), (X(centloc(1), centloc(2)) +lof2), (X(centloc(1), centloc(2)) -lof2)];
                Ypi     = [Y(centloc(1), centloc(2)), (Y(centloc(1),...
                    centloc(2)) -lof2), (Y(centloc(1), centloc(2)) +lof2), (Y(centloc(1), centloc(2)) -lof1), (Y(centloc(1), centloc(2)) +lof1)];
                Xp      = Xpi;
                Yp      = Ypi;

                % If you want to savesome information
                %lof1save(mm,:) = lof1;
                %lof2save(mm,:) = lof2;
                %xpisave(mm,:)  = Xpi;
                %ypisave(mm,:)  = Ypi;

                for nn = 1:nt

                    % Integrating dx/dt=u,dy/dt=v over a time step dt using Improved Euler Integration
                    % up, vp - velocities at positions at beginning of time-step
                    up=interp2(Y,X,u_interp,Yp,Xp);
                    vp=interp2(Y,X,v_interp,Yp,Xp);

                    % Find NaN values in up/vp during advection
%                     tmp_up(:,jj,nn) = interp2(Y,X,u,Yp,Xp);
%                     tmp_vp(:,jj,nn) = interp2(Y,X,v,Yp,Xp);
%                     tmp_xp(:,jj,nn) = Xp;
%                     tmp_yp(:,jj,nn) = Yp;

                    % First step of Improved Euler is Forward Euler (this is the final position of the velocity data for the entire time window)
                    Xpstar=Xp+up*dt;
                    Ypstar=Yp+vp*dt;

                    % Second step of Improved Euler
                    upstar=interp2(Y,X,u_interp,Ypstar,Xpstar);
                    vpstar=interp2(Y,X,v_interp,Ypstar,Xpstar);

                    % Final positions after each Improved Euler Time Step
                    Xp=Xp+(up+upstar)/2*dt;
                    Yp=Yp+(vp+vpstar)/2*dt;

                   
                    %sprintf('timestep %d out of %d, Moving on to the next timestep',nn, length(nt))

                end %nn

                % Save Xp/Yp values for each rotation

%                 alltmp_up(mm,:,:,:) = tmp_up;
%                 alltmp_vp(mm,:,:,:) = tmp_vp;
%                 alltmp_xp(mm,:,:,:) = tmp_xp;
%                 alltmp_yp(mm,:,:,:) = tmp_yp;

                % Final lengths based on final coordinates
                a1f        = sqrt((Xp(5)-Xp(1))^2+(Yp(5)-Yp(1))^2);
                a2f        = sqrt((Xp(4)-Xp(1))^2+(Yp(4)-Yp(1))^2);
                b1f        = sqrt((Xp(5)-Xp(2))^2+(Yp(5)-Yp(2))^2);
                b2f        = sqrt((Xp(4)-Xp(3))^2+(Yp(4)-Yp(3))^2);
                c1f        = sqrt((Xp(2)-Xp(1))^2+(Yp(2)-Yp(1))^2);
                c2f        = sqrt((Xp(3)-Xp(1))^2+(Yp(3)-Yp(1))^2);
                d1f        = sqrt((Xp(4)-Xp(2))^2+(Yp(4)-Yp(2))^2);
                d2f        = sqrt((Xp(5)-Xp(3))^2+(Yp(5)-Yp(3))^2);

                % Calculating logarithmic strain rates
                tott       = nt*dt; % Total time for virtual advection
                edot0      = 0.5/tott*(log(a1f/a10)+log(a2f/a10));
                edot45     = 0.5/tott*(log(b1f/b10)+log(b2f/b10));
                edot90     = 0.5/tott*(log(c1f/c10)+log(c2f/c10));
                edot135    = 0.5/tott*(log(d1f/d10)+log(d2f/d10));

                % a legs
                ca         = cosd(alpha);sa=sind(alpha);
                % b legs
                cam45      = cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%cos(alpha-45)
                sam45      = -cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%sin(alpha-45)
                % c legs
                cam90      = sind(alpha);%cos(alpha-90)
                sam90      = -cosd(alpha);%sin(alpha-90)
                % d legs
                cap45      = cosd(alpha)/sqrt(2)-sind(alpha)/sqrt(2);%cos(alpha+45)
                sap45      = cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%sin(alpha+45)
                %
                A          = [ca^2,2*ca*sa,sa^2;cam45^2,2*cam45*sam45,sam45^2;...
                    cam90^2,2*cam90*sam90,sam90^2;cap45^2,2*cap45*sap45,sap45^2];
                lsqmatrix  = inv(A'*A)*A';
                edotvector = lsqmatrix*[edot0;edot45;edot90;edot135];

                % Components of the strain rate tensor
%                 exxdot     = edotvector(1);
%                 exydot     = edotvector(2);
%                 eyydot     = edotvector(3);
                
                xpisave(:,mm)       = Xpi;
                ypisave(:,mm)       = Ypi;
                xpsave(:,mm)        = Xp;
                ypsave(:,mm)        = Yp;
                eresults(:,mm)      =...
                                        [edotvector(1),edotvector(2),edotvector(3)]; % Write estimated components into row of eresults
                %clear Yp Xp Xpi Ypi edotvector
            end %mm each alphas
                ls(jj).eresults(:,:,kk) = eresults;
                ls(jj).xpi(:,:,kk) = xpisave;
                ls(jj).ypi(:,:,kk)  = ypisave;
                ls(jj).xp(:,:,kk) = xpsave;
                ls(jj).yp(:,:,kk) = ypsave;
                ls(jj).lengthscale(kk) = couplinglengthscale(kk);

                clear xpisave ypisave xpsave ypsave eresults
            %sprintf('Finished with length scale %d out of %d, Moving on to the next lengthscale',kk, length(couplinglengthscale))
        end %kk lengthscales 
            
            
    end %jj velocity perturbations
            strainrates(ii).vp      = ls;
            strainrates(ii).crossoverROI   = is2roi(ii,:);
            strainrates(ii).crossoverthick = roithick(ii);
   %           strainrates(ii).final_xoverthick =
            strainrates(ii).roibederror = roibederror(ii);

    progress_xover = sprintf('Finished with crossover %d out of %d, Moving on to the next crossover',ii, length(is2roi));
    disp(progress_xover)
    toc
end %ii crossover

% Using the centered is2 location, generate a list (XY_list) of x and y coordinates matching the five points of the strain diamond coordinates (Xpi,Ypi) and the final positions (Xp,Yp) for each rotation and length scale.
% XY_list = struct('xpi',[],'ypi',[],'xp',[],'yp',[]);
% for ii = 1:size(is2roi,1) %was centered_is2_locations
%     for jj = 1:size(couplinglengthscale,2)
%         for kk = 1:size(numsim,2) % numsim was velperturb
%             XY_list(ii).xpi(:,:,jj,kk) = ls(jj).xpi(:,:,kk);
%             XY_list(ii).ypi(:,:,jj,kk) = ls(jj).ypi(:,:,kk);
%             XY_list(ii).xp(:,:,jj,kk) = ls(jj).xp(:,:,kk);
%             XY_list(ii).yp(:,:,jj,kk) = ls(jj).yp(:,:,kk);
%             strainrates(ii).vp(jj,kk).eresults = ls(jj).eresults(:,:,kk);
%             strainrates(ii).vp(jj,kk).lengthscale = ls(jj).lengthscale(kk);
%         end
%     end
% end
% Do this for each crossover. Turn those X and Y coordinates into a list of strain rates for each location point. Save the strain rates, the coordinates, and the length scales for each crossover. 
% Turn the X and Y coordinate values in XY_list into latitude and longitude values and save those as well.



output.strainrates                 = strainrates;
output.velocityinfo.pixelsize      = vel.pixelsize;
output.velocityinfo.vel_starttime  = vel.start_time;
output.velocityinfo.vel_endtime    = vel.end_time;
output.velocityinfo.e_vel          = vel.e_vel;
output.velocityinfo.n_vel          = vel.n_vel;
output.velocityinfo.x_vel          = vel.x;
output.velocityinfo.y_vel          = vel.y;
% output(i).thickinfo.thickness         = thick.thickness;
% output(i).thickinfo.thickerror        = thick.errbed;
% output(i).thickinfo.thickx            = thick.x;
% output(i).thickinfo.thicky            = thick.y;
% output(i).thickinfo.source            = thick.source;

disp('Finished with velocity scene, outputs ready')


if saveoutputs == 1
    fname = sprintf('SR_for_%s_Cycle_%d_xover_%6.1f_%6.1f.mat',Region,Cycle,is2roi(ii,:)); %is2roi was centered_is2_locations(ii,:)
    save(fname, 'strainrates')
end












