function [roithick, roibederror, velrow, velcol, vel] = is2roiinfo(is2roi, thick, vel, largeboundingboxsize)

xmax = max(is2roi(:,1), [], "all") + largeboundingboxsize*1000; % do we want a smaller box still?
xmin = min(is2roi(:,1), [], "all") - largeboundingboxsize*1000;
ymax = max(is2roi(:,2), [], "all") + largeboundingboxsize*1000; 
ymin = min(is2roi(:,2), [], "all") - largeboundingboxsize*1000;

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
tmpthickxy(:,4) = thick.errbed(:);
tmpthickxy(any(isnan(tmpthickxy), 2), :) = [];

for ii = 1:length(is2roi)
    [index2(ii), distance2(ii)] = knnsearch(tmpthickxy(:,1:2), is2roi(ii,:), 'K', 1,'distance', 'euclidean');
end

roithick = tmpthickxy(index2,3); %in METERS
roibederror = tmpthickxy(index2,4); %in METERS 

%%%%%%%%

%for kk = 1:length(vel)
%     vel.x_tmp = cell2mat(vel(1,1).x); %%%%changed to x_tmp bc cell/double conversion
%     vel.y_tmp = cell2mat(vel(1,1).y); %%%%changed to y_tmp bc cell/double conversion
%     tmpvelx = double(vel.x_tmp); % for cells
%     tmpvely = double(vel.y_tmp); % for cells

for i = 1:length(vel)
    vel_struct = vel{i, 1}; % Extract the nested structure - Change {1,1} to {i,1} for multiple scenes
    
    % Clear previous iteration's data
    %clear x_values y_values
    
    x_values = vel_struct.x; % Access the 'x' field
    y_values = vel_struct.y; % Access the 'y' field
    
    tmpvelx = double(x_values);
    tmpvely = double(y_values);

    clear vel_struct x_values y_values
%     tmpvelx = double(vel.x); % remove if using cells
%     tmpvely = double(vel.y); % remove if using cells
    tmpvelx((tmpvelx<xmin)) = NaN;
    tmpvelx((tmpvelx>xmax)) = NaN;
    tmpvely((tmpvely<ymin)) = NaN;
    tmpvely((tmpvely>ymax)) = NaN;

    tmponev = tmpvelx + tmpvely;
    tmponev(~isnan(tmponev)) = 1;

    indicies1 = find(tmponev(:)==1);
    [upleft(1,1), upleft(1,2)] = ind2sub(size(tmponev),indicies1(1)) ;
    [botright(1,1), botright(1,2)] = ind2sub(size(tmponev),indicies1(end));
    
    vel{i, 1}.y = vel{i, 1}.y(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.x = vel{i, 1}.x(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.n_vel = vel{i, 1}.n_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.e_vel = vel{i, 1}.e_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.n_vel_std = vel{i, 1}.n_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.e_vel_std = vel{i, 1}.e_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.tot_vel = vel{i, 1}.tot_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.tot_vel_std = vel{i, 1}.tot_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.n_vel_base = vel{i, 1}.n_vel_base(upleft(1):botright(1), upleft(2): botright(2)) ;
    vel{i, 1}.e_vel_base = vel{i, 1}.e_vel_base(upleft(1):botright(1), upleft(2): botright(2)) ;
    

    tmpvelxy(:,1) = vel{i,1}.x(:);
    tmpvelxy(:,2) = vel{i,1}.y(:);

    for ii = 1:length(is2roi)
        [index2(ii), distance2(ii)] = knnsearch(tmpvelxy(:,1:2), is2roi(ii,:), 'K', 1,'distance', 'euclidean');
        [velrow(ii), velcol(ii)] = ind2sub(size(vel{i,1}.x),index2(ii));
    end

    velrow = velrow';
    velcol = velcol';
    vel

    %     vel.y = vel.y(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.x = vel.x(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.n_vel = vel.n_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.e_vel = vel.e_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.n_vel_std = vel.n_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.e_vel_std = vel.e_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.tot_vel = vel.tot_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.tot_vel_std = vel.tot_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.n_vel_base = vel.n_vel_base(upleft(1):botright(1), upleft(2): botright(2)) ;
%     vel.e_vel_base = vel.e_vel_base(upleft(1):botright(1), upleft(2): botright(2)) ;
%     %vel.vert_vel = vel.vert_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
%     % For uncertainty, incorporate e_vel_std and n_vel_std here from readvelocitydata

end
