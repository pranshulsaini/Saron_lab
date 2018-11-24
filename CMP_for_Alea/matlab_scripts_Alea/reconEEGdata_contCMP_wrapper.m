% Last modified before Pranshul: 24th Feb 2015

% This code is written as a wrapper for reconEEGdata_forCMPdata.m
% Input to this function is an .xlsx file containing one column that lists
% the name of the subject whose data is to be reconstructed. 
%%
% Expectations:
%(1) Subject has been run through SOBI and sources have been voted on
%(2) Subject's original preprocessed data (.dat, .evt, .sfp, and .generic file) is in C:\forSOBI_RestingData\
%(3) Source votes for each subject are in .xlsx format named
%'RST#####_voteforrecon', which is saved in the same directory as the SMART
%output. The excel file should follow the format described in the protocol,
%with column A containing source numbers and column B containing a vote of
%"norm" for presumed neural sources and "emg" for presumed noise sources.

%==================================================================
function reconEEGdata_contCMP_wrapper(xlsfilename)


    % adding path
    addpath('C:\Alea\matlab_scripts');
    addpath('C:\Alea\matlab_scripts\extra_required_codes');

    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
    
    % r =   {' FILES TO BE RECONSTRUCTED'}       %(P)
    %       {'STR00412'                  }

    
    % start from i = 2, since the first row is titles
    for i = 2:1:size(r,1)        % size(r,1) is (one + the number of subjects' data to be reconstructed )  (P)
        nameoffile = r{i,1};     % 'STR00412'  (P)
        outloc = strcat('C:\Alea\CMP_SOBI_Output\',nameoffile,'_forSOBI\');
        eegfile = strcat([nameoffile, '_forSOBI']); %  'STR00412_aux_NoBad_AvgRef_forSOBI'. No extension (P)        
        eegdir = strrep(eegfile, '.dat','');  % No effect. Stays 'STR00412_aux_NoBad_AvgRef_forSOBI'. Probably that's what we want (P)
        EOdataloc = strcat(outloc, [eegdir, '_Sources_CMP\']);  % 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\STR00412_aux_NoBad_AvgRef_forSOBI_Sources_STR\'  (P)
        SmartOutloc = strcat(EOdataloc, [eegfile,'smarter_CMP\']); %   ''C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\STR00412_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR00412_aux_NoBad_AvgRef_forSOBIsmarter_STR'
        createOutFromXLS([outloc, eegfile, '_voteforrecon.xlsx'],[EOdataloc,[eegdir, '.out']]);  % takes an excel file (of votes) as input and creates a similar .out file for reconstructing the EEG data (P)
        % The directory of out file would be: 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\STR00412_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR00412_aux_NoBad_AvgRef_forSOBI.out'
        
        outfile = [eegdir, '.out']; % 'STR00412_aux_NoBad_AvgRef_forSOBI.out' This is the file created by createOutFromXLS (P)
%         CreateNoiseReconFile([SmartOutloc, nameoffile, '_voteforrecon_EO.xlsx'],[EOdataloc,[eegdir, '_EyesOpenNoise.out']]);
%         noiseoutfile = [eegdir, '_EyesOpenNoise.out'];
%         %write the reconstructed data
        reconEEGdata_forCMPdata(eegdir, 1, EOdataloc,  outfile); % 
        % writes clean reconstructed data in terms of .dat, .generic and .evt files

        reconEEGdata_forCMPdata(eegdir, 2, EOdataloc,  outfile); % 
        % writes noise data (also called reconInverse) in terms of .dat, .generic and .evt files

    end
end


