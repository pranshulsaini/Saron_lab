% This code is written as a wrapper for tSOBI_CPT_{etc}.m This file will
% process the list of filenames in the xlsfilename submitted to the
% function command.


function GenerateSourceProjections_wrapper(xlsfilename) 

    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
    
    % location of DAT, EVT and GENERIC files
    datloc = 'H:\CPT_SOBI_data\';
    
    % location of sfp files
%    sfploc = 'H:\CPT_SOBI_data\';
    
    % start from i = 2, since the first row is header
    for i = 2:1:size(r,1)
        eegfilename = strrep(r{i,1},'.dat','');
        eegfile = strcat(datloc, r{i,1});
%        evtfile = strcat(datloc, strrep(r{i,1}, '.dat', '.evt'));
%        sfpfile = strcat(sfploc, strrep(r{i,2}, '.dat', '.evt'));

        % create folder for each subject
        output = strrep(r{i,1}, '.dat', ''); 
        outputDir = strcat('H:\CPT_SOBI_Output\',output);
        if exist(outputDir,'dir') ~= 7
           fprintf(1, '\nDirectory for the subject (%s) doesn''t exist, Creating one...',output);
           mkdir(outputDir);
        else 
           fprintf(1, '\nArtifacts dir already exists, Overwriting...'); 
        end
        cd(outputDir);          
        chunks = 1;
            RegenerateSourceProjections(eegfilename, eegfile,outputDir);    %generate source projections
       
    end
end