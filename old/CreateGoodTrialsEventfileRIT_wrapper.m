% This code is written as a wrapper for CreateGoodTrialsEventfile.m
% Input to this function is an .xlsx file (GoodTs.xlsx)containing one column that lists
% the name of the subject whose data is to be reconstructed. (Example:
% CPT00321)

function CreateGoodTrialsEventfile_wrapper(xlsfilename)

    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
     
       % start from i = 2, since the first row is titles
    for i = 2:1:size(r,1)
        nameoffile = r{i,1};
       CreateGoodTrialsEventfileRIT(nameoffile);
   %    disp(['Good Trials event file generated for ', nameoffile]);
    end
end
