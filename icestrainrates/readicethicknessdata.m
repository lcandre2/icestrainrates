function t = readicethicknessdata(datasource, path1, filetype)

if datasource == 1
    
    location1 = [path1, filetype];
    
    files        = dir(location1); %will load all the nc4 files in the folder
    files        = struct2cell(files);
    files        = files(1, 1:end);
    files        = files';

   files2 = [path1, files{1}];
    %This information is specific to promice data
    t.x      = ncread(files2, 'x');
    t.y      = ncread(files2, 'y');
    [t.x, t.y] = meshgrid(t.x, t.y);
    t.x = t.x'; t.y=t.y';

    for ii = 1:length(files)
        clear files2
        files2 = [path1, files{ii}];
        t.mask(:,:,ii)  = ncread(files2, 'mask');
        t.thickness(:,:,ii)  = ncread(files2, 'thickness');

    end


elseif datasource == 2
%

elseif datasource == 3



else 
    disp('The code does not understand your choice of ice thickness data: Abort!')

end 