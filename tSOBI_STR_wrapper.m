% This code is written as a wrapper for tSOBI_STR_{etc}.m This file will
% process the list of filenames in the xlsfilename submitted to the
% function command. It has been modified by Pranshul Saini

% Last modified before Pranshul: 2nd April 2015

function tSOBI_STR_wrapper(xlsfilename)
    time_beg = cputime;

    addpath('C:\Users\plsaini\Box Sync\Stroop\Temp\matlab_scripts');
    addpath('C:\Users\plsaini\Box Sync\Stroop\Temp\matlab_scripts\extra_required_codes');
    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
    
    % location of DAT, EVT and GENERIC files
    datloc = 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR_SOBI_data\';
    
    % location of sfp files
    sfploc = 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR_SOBI_data\';
    
    % start from i = 2, since the first row is header
    for i = 2:1:size(r,1)
        eegfilename = strrep(r{i,1},'.dat','');  % replaces 'dat' with nothing. So, eegfilename = 'STR00412_aux_NoBad_AvgRef_forSOBI'
        eegfile = strcat(datloc, r{i,1});  % dat file
        evtfile = strcat(datloc, strrep(r{i,1}, 'aux_NoBad_AvgRef_forSOBI.dat', '1stPass.evt'));  % replaces aux_NoBad_AvgRef_forSOBI.dat with 1stPass.evt
        %evtfile = strrep(r{i,1}, '.dat', '.evt');
        sfpfile = strcat(sfploc, strrep(r{i,2}, '.dat', '.sfp')); % replaces.dat with .sfp

        % create folder for each subject
        output = strrep(r{i,1}, '.dat', ''); % name of the folder
        outputDir = strcat('C:\Users\plsaini\Box Sync\Stroop\Temp\',output);
        if exist(outputDir,'dir') ~= 7   % it returns 7 if outputDir exists as a folder
           fprintf(1, '\nA directory for the subject (%s) doesn''t exist, Creating one...',output);
           mkdir(outputDir);
        else 
           fprintf(1, '\nA directory for the subject already exists, Overwriting...'); 
        end
        cd(outputDir); % enters this folder      
        chunks = 1;  % I do not what it means but I know it remains 1 or CPT, as well as Stroop
%       tSOBI_CPT_And_Single_Source_Recon_2700ms_difout(eegfilename, eegfile, evtfile, sfpfile, outputDir, chunks);    %run SOBI
     
        tSOBI_STR_betterSMART(eegfilename, eegfile, evtfile, sfpfile, outputDir, chunks);    %run SOBI
       
    end
    
    e = cputime - time_beg;
    fprintf('The time taken for the simulation is %f', e);
end