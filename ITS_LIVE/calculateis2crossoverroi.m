function is2roi = calculateis2crossoverroi(pathname, filetype)

   location1     = [pathname, filetype];
    files        = dir(location1); %will load all the nc4 files in the folder
    files        = struct2cell(files);
    files        = files(1, 1:end);
    files        = files';

    files2 = [pathname, files{1}];

laser_xy = [NaN, NaN];
% Find the indicies of the points in a cluster
for ii = 1:length(files)
     tmp1  = textread(files2);
     laser_xy = [laser_xy; tmp1];
end

laser_xy = laser_xy(2:end, :);

for ii = 1:length(laser_xy)
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
        is2roi(ii,jj) = mean(laser_xy(tmp,jj),1, 'omitnan');
    end
end

