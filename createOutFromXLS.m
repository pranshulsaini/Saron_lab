% Last modified before Pranshul: 11th April 2012

% This code is written for taking an excel file (of votes) as input and
% creating an .out file for reconstructing the EEG data.
%
% USAGE: createOutFromXLS(xlsfilename, nameOfOutFileWithPath)
%  e.g. createOutFromXLS('/Users/Shared/autism/votes_example.xls',
%  '/Users/Shared/autism/votes_example.out')
%  This function call will create a new .out file named votes_example in
%  the folder /Users/Shared/autism
% 
% Written by Manish Saggar (mishu@cs.utexas.edu) on 2/15/2011

function createOutFromXLS(xlsfile, output)
    [n,t,r] = xlsread(xlsfile);
    fid = fopen(output,'w');
    if fid == -1
        fprintf(2, 'Error creating out file %s, make sure path exists.', output);
        return;
    end
    
    for i = 1:1:size(r,1)
        % check if there is a header            In that case, it skips and goes to the next row (P)
        if ~isnumeric(r{i,1})         
            continue;
        end
   
        fprintf(fid, '%d %s\n', r{i,1}, r{i,2});
    end
    fclose(fid);
end