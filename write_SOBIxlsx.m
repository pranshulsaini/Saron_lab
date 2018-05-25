function write_SOBIxlsx(inloc, outloc, name)

    runsDir = dir(inloc); %# Get the data for the current directory
    dirIndex = [runsDir.isdir];  %# Find the index for directories
    fileList = {runsDir(~dirIndex).name}';

    filename = strcat(outloc, name);

    col_header={'DAT FILES','SFP FILES'};  


    x = strings(length(fileList),1);
    for i = 1:length(fileList)
        x(i) = fileList{i};
    end
    
    data = strings(length(fileList),2);


    tmp = strrep(x,'.bdf','write');
    data(:,1) = tmp + '_aux_NoBad_AvgRef_forSOBI.dat';
    data(:,2) = tmp + '_aux_NoBad_AvgRef_forSOBI.sfp';


    xlswrite(filename,col_header,'Sheet1','A1');     %Write column header

    xlswrite(filename,data,'Sheet1','A2');   % write data

end



