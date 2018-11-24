% This code is written as a wrapper for CreateUniqueIDs.m
% Input to this function is an .xlsx file (UIDS.xlsx)containing one column that lists
% the name of the subject whose data is to be reconstructed. (Example:
% CPT00321)

function CreateUniqueIDsRIT_wrapper(xlsfilename)

    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
     
       % start from i = 2, since the first row is titles
    for i = 2:1:size(r,1)
        nameoffile = r{i,1};
       CreateUniqueIDsRIT(nameoffile);
       disp(['Unique IDs generated for ', nameoffile]);
    end
end
