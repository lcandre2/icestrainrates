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
vel.n_vel_std = vel.n_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.e_vel_std = vel.e_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.tot_vel = vel.tot_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
% vel.vert_vel = vel.vert_vel(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.tot_vel_std = vel.tot_vel_std(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.n_vel_base = vel.n_vel_base(upleft(1):botright(1), upleft(2): botright(2)) ;
vel.e_vel_base = vel.e_vel_base(upleft(1):botright(1), upleft(2): botright(2)) ;
% For uncertainty, incorporate e_vel_std and n_vel_std here from readvelocitydata

tmpvelxy(:,1) = vel.x (:); 
tmpvelxy(:,2) = vel.y(:);

for ii = 1:length(is2roi)
    [index2(ii), distance2(ii)] = knnsearch(tmpvelxy(:,1:2), is2roi(ii,:), 'K', 1,'distance', 'euclidean');
    [velrow(ii), velcol(ii)] = ind2sub(size(vel.x),index2(ii)); 
end

velrow = velrow';
velcol = velcol';
vel
