% This code is written as a wrapper for tSOBI_CPT_{etc}.m This file will
% process the list of filenames in the xlsfilename submitted to the
% function command.


function tSOBI_CPT_wrapper(xlsfilename) 

    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
    
    % location of DAT, EVT and GENERIC files
    datloc = '\\Dss02544-cmb-d\h\CPT_SOBI_data\';
    
    % location of sfp files
    sfploc = '\\Dss02544-cmb-d\h\CPT_SOBI_data\';
    
    % start from i = 2, since the first row is header
    for i = 2:1:size(r,1)
        eegfilename = strrep(r{i,1},'.dat','');  % replaces 'dat' with nothing. So, eegfilename = 'CPT03811_forSOBI'
        eegfile = strcat(datloc, r{i,1});  % dat file
        evtfile = strcat(datloc, strrep(r{i,1}, 'forSOBI.dat', 'GoodTrials.evt'));  % replaces forSOBI.dat with GoodTrials.evt
        %evtfile = strrep(r{i,1}, '.dat', '.evt');
        sfpfile = strcat(sfploc, strrep(r{i,2}, '.dat', '.evt')); % replaces.dat with .evt

        % create folder for each subject
        output = strrep(r{i,1}, '.dat', ''); % name of the folder
        outputDir = strcat('\\Dss02544-cmb-d\h\CPT_SOBI_Output\',output);
        if exist(outputDir,'dir') ~= 7   % it returns 7 if outputDir exists as a folder
           fprintf(1, '\nA directory for the subject (%s) doesn''t exist, Creating one...',output);
           mkdir(outputDir);
        else 
           fprintf(1, '\nA directory for the subject already exists, Overwriting...'); 
        end
        cd(outputDir); % enters this folder      
        chunks = 1;  % I do not what it means
%       tSOBI_CPT_And_Single_Source_Recon_2700ms_difout(eegfilename, eegfile, evtfile, sfpfile, outputDir, chunks);    %run SOBI
        tSOBI_CPT_betterSMART(eegfilename, eegfile, evtfile, sfpfile, outputDir, chunks);    %run SOBI
       
    end
end