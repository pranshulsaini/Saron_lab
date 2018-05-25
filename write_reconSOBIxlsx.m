function read_filenames(inloc, outloc, name)

    runsDir = dir(inloc); %# Get the data for the current directory
    dirIndex = [runsDir.isdir];  %# Find the index for directories
    fileList = {runsDir(~dirIndex).name}';

    filename = strcat(outloc, name);

    col_header={'FILES TO BE RECONSTRUCTED'};  

    x = strings(length(fileList),1);

    for i = 1:length(fileList)
        x(i) = fileList{i};
    end

    data = strrep(x,'.bdf','');


    xlswrite(filename,col_header,'Sheet1','A1');     %Write column header

    xlswrite(filename,data,'Sheet1','A2');   % write data


end


