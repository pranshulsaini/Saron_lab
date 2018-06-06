function read_filenames(inloc, outloc, name)

% remember to add '\' in the end of inloc and outloc

runsDir = dir(inloc); %# Get the data for the current directory
dirIndex = [runsDir.isdir];  %# Find the index for directories
fileList = {runsDir(~dirIndex).name}';

filename = strcat(outloc, name);
fid = fopen(filename,'w');
for i = 1:length(fileList)
  fprintf(fid,'%s\n',fileList{i});
end

fclose(fid);


end

