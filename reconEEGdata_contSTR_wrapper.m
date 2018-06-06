% Last modified before Pranshul: 24th Feb 2015

% This code is written as a wrapper for recondEEGdata_forRestdata.m
% Input to this function is an .xlsx file containing one column that lists
% the name of the subject whose data is to be reconstructed. (Example:
% RST00321)
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
function reconEEGdata_contSTR_wrapper(xlsfilename)


    % adding path
    addpath('C:\Users\plsaini\Box Sync\Stroop\matlab_scripts');
    addpath('C:\Users\plsaini\Box Sync\Stroop\extra_matlab');

    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
    
    % r =   {' FILES TO BE RECONSTRUCTED'}       %(P)
    %       {'STR00412'                  }

     
    % location of DAT, EVT and GENERIC files
    
    %outloc = 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\';
    
    % start from i = 2, since the first row is titles
    for i = 2:1:size(r,1)        % size(r,1) is (one + the number of subjects' data to be reconstructed )  (P)
        nameoffile = r{i,1};     % 'STR00412'  (P)
        outloc = strcat('C:\Users\plsaini\Box Sync\Stroop\STR_SOBI_Output\',nameoffile,'_aux_NoBad_AvgRef_forSOBI\');
        eegfile = strcat([nameoffile, '_aux_NoBad_AvgRef_forSOBI']); %  'STR00412_aux_NoBad_AvgRef_forSOBI'. No extension (P)        
        evtfile = strrep(eegfile, '.dat', '.evt');  % Does not work. Stays 'STR00412_aux_NoBad_AvgRef_forSOBI'. No extension added. I should comment it out  (P)
        sfpfile = strrep(eegfile, '.dat', '.sfp');  % Does not work. Stays 'STR00412_aux_NoBad_AvgRef_forSOBI'. No extension added. I should comment it out  (P)
        eegdir = strrep(eegfile, '.dat','');  % No effect. Stays 'STR00412_aux_NoBad_AvgRef_forSOBI'. Probably that's what we want (P)
        EOdataloc = strcat(outloc, [eegdir, '_Sources_STR\']);  % 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\STR00412_aux_NoBad_AvgRef_forSOBI_Sources_STR\'  (P)
        SmartOutloc = strcat(EOdataloc, [eegfile,'smarter_STR\']); %   ''C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\STR00412_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR00412_aux_NoBad_AvgRef_forSOBIsmarter_STR'
        createOutFromXLS([outloc, eegfile, '_voteforrecon.xlsx'],[EOdataloc,[eegdir, '.out']]);  % takes an excel file (of votes) as input and creates a similar .out file for reconstructing the EEG data (P)
        % The directory of out file would be: 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\STR00412_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR00412_aux_NoBad_AvgRef_forSOBI.out'
        
        outfile = [eegdir, '.out']; % 'STR00412_aux_NoBad_AvgRef_forSOBI.out' This is the file created by createOutFromXLS (P)
%         CreateNoiseReconFile([SmartOutloc, nameoffile, '_voteforrecon_EO.xlsx'],[EOdataloc,[eegdir, '_EyesOpenNoise.out']]);
%         noiseoutfile = [eegdir, '_EyesOpenNoise.out'];
%         %write the reconstructed data
        reconEEGdata_forSTRdata(eegdir, 1, EOdataloc,  outfile); % 
        % writes clean reconstructed data in terms of .dat, .generic and .evt files

        reconEEGdata_forSTRdata(eegdir, 2, EOdataloc,  outfile); % 
        % writes noise data (also called reconInverse) in terms of .dat, .generic and .evt files


        %  reconNoise_forRestdata(eegdir, 1, EOdataloc, outfile, noiseoutfile);
        % reconNoise_forRestdata(eegdir, 2, EOdataloc, outfile, noiseoutfile);
        IntegrateResponseTrigs(eegfile);
    end
end


function IntegrateResponseTrigs(nameoffile)  % nameoffile = 'STR00412_aux_NoBad_AvgRef_forSOBI'

    reconevtfile = strcat('recon_',nameoffile, '.evt'); %original              This is the event file which got created in the file reconEEGdata_forSTRdata.m (P)
    %reconevtfile = strcat('recon_',nameoffile, '_CPT.evt'); %changed for the test data. change back to orig for normal files
    %eegfile = strcat(nameoffile, '_epoched_forSOBI'); %original
    eegfile = nameoffile; %changed for the test data. change back to orig for normal files
    resptrigevtfile = strcat(eegfile(1:8), '_1stPass.evt');    % removed 'for_SOBI' and added '_GoodTrials.evt' (P)
    %outloc = 'C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\';   
    outloc = strcat('C:\Users\plsaini\Box Sync\Stroop\STR_SOBI_Output\',nameoffile,'\');
    recondataloc = strcat(outloc, [nameoffile, '_Sources_STR\'], [nameoffile, '_recon\']); % folder for reconstructed data (P)  'C:\Users\plsaini\Box Sync\Stroop\Temp\STR00412_aux_NoBad_AvgRef_forSOBI\STR00412_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR00412_aux_NoBad_AvgRef_forSOBI_recon\'; 
    Trigdataloc = 'C:\Users\plsaini\Box Sync\Stroop\STR_SOBI_data\';
 
    
    %This is the preSOBI GoodTrial event file which has the information about the original response triggers timings (P)
    filename = strcat(Trigdataloc,resptrigevtfile);
    fid = fopen(filename,'r');
    origevt = textscan(fid,'%d %d %d %s %d %s %s','headerlines',1);   % will scan the whole text of the event file. First line will be assumed as heading
    fclose(fid);
   
    n_rows_orig = size(origevt{1,1},1);
    origtrig = origevt{1,3}(:);  % original triggers from 1stPass event file
    origtmu = origevt{1,1}(:);
    
    filename = strcat(recondataloc,reconevtfile);
    fid = fopen(filename,'r');
    currevt = textscan(fid,'%d %d %d %s','headerlines',1);   % will scan the whole text of the event file. First line will be assumed as heading
    fclose(fid);
    
    currtmu = currevt{1,1}(:);  % current timestamp values
    currtrig = currevt{1,3}(:);  % current trigger values
    currcomm = currevt{1,4}(:);  % current comment values
    
    n_rows_recon = 2*size(currevt{1,1},1);   %  2 is used because we want to accommodate response triggers as well
    recontmu = zeros(n_rows_recon,1);   % they will store final timestamps for the reconstructed event file
    recontrig = zeros(n_rows_recon,1);   % they will store final triggers for the reconstructed event file
    reconcode = ones(n_rows_recon,1);  % the second column
    
    reconcomm =  strings(n_rows_recon,1);   % string is important to write to the event file
    reconcomm(1:2:n_rows_recon) =  currcomm;  % comments will be the trial number
    reconcomm(2:2:n_rows_recon) =  currcomm;  % comments will be the trial number
 
    recontmu(1:2:end) = currtmu + .25e6;
    recontrig(1:2:end) = currtrig;
    
    j = 2;
    for i = 1:n_rows_orig    % 
        if (origtrig(i) > 1000) % This will let pass only the stimulus triggers
            if (i < n_rows_orig)   % just to be safe for end trial
                RT = origtmu(i+1) - origtmu(i);      % assuming the next trigger will be a response trigger. This will be checked in the next line
                if (origtrig(i+1) < 1000) && (RT< 1750000)   % the second condition makes sure that they are consequent
                    recontmu(j) = recontmu(j-1) + RT;    % so this is the time stamp for the response
                    recontrig(j) = origtrig(i+1);  % just the corresponding tiggers
                else
                    recontmu(j) = recontmu(j-1) + 1500000;   % added as an arbitrary number
                    recontrig(j) = 500;   % labeling as miss
                end
            else
                recontmu(j) = recontmu(j-1) + 1500000;
                recontrig(j) = 500;   % labelling as miss
            end
        j = j + 2;    
        end    
    end
            
    
    % writing the integrated event file
    filename =strcat(recondataloc,'recon_',eegfile,'_withRespTrigs.evt');
    fid = fopen(filename,'w');
    
    if fid == -1
        fprintf(1,'Error creating event file with integrated response triggers.\n');
    end
    
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    
    for j = 1:n_rows_recon
        fprintf(fid, '%d\t%d\t%d\t%s\n', recontmu(j), reconcode(j), recontrig(j), reconcomm(j));
    end
    
    fclose(fid);
    
    disp('Response Triggers successfully integrated.');
    
end              



function CreateNoiseReconFile(xlsfilename, output)
%read the Excel file
[n,t,r] = xlsread(xlsfilename);

%creat new file with channels marked as norm changed to emg and vice versa
fid = fopen(output,'w');
if fid == -1
        fprintf(2, 'Error creating out file %s, make sure path exists.', output);
        return;
end
    for i = 1:1:size(r,1)
        if strcmp(r(i,2),'norm') == 1
            correction = 'emg';
        else
            correction = 'norm';
        end
            fprintf(fid, '%d %s\n', r{i,1}, correction);
    end
    fclose(fid)
end