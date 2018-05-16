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
function reconEEGdata_contCPT_wrapper(xlsfilename)

    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
    
    % r =   {' FILES TO BE RECONSTRUCTED'}       %(P)
    %       {'CPT03811'                  }

    
     
    % location of DAT, EVT and GENERIC files
    
    
    %OGdatloc = '\\Dss02721-cmb-d\c\forSOBI_RestingData\'; %location of original data files sent through SOBI
    outloc = '\\Dss02544-cmb-d\h\CPT_SOBI_Output\';
    % start from i = 2, since the first row is titles
    for i = 2:1:size(r,1)        % size(r,1) is (one + the number of subjects' data to be reconstructed )  (P)
        nameoffile = r{i,1};     % 'CPT03811'  (P)
        eegfile = strcat([nameoffile, '_forSOBI']); %  'CPT03811_forSOBI'. No extension (P)        
        evtfile = strrep(eegfile, '.dat', '.evt');  % Does not work. Stays 'CPT03811_forSOBI'. No extension added. I should comment it out  (P)
        sfpfile = strrep(eegfile, '.dat', '.sfp');  % Does not work. Stays 'CPT03811_forSOBI'. No extension added. I should comment it out  (P)
        eegdir = strrep(eegfile, '.dat','');  % No effect. Stays 'CPT03811_forSOBI'. Probably that's what we want (P)
        EOdataloc = strcat(outloc, eegdir,'\', [eegdir, '_Sources_CPT\']);  % '\\Dss02544-cmb-d\h\CPT_SOBI_Output\CPT03811_forSOBI\CPT03811_forSOBI_Sources_CPT\'  (P)
        SmartOutloc = strcat(EOdataloc, [eegfile,'smarter_CPT\']); %   '\\Dss02544-cmb-d\h\CPT_SOBI_Output\CPT03811_forSOBI\CPT03811_forSOBI_Sources_CPT\CPT03811_forSOBIsmarter_CPT\'
        createOutFromXLS([SmartOutloc, eegfile, '_voteforrecon.xlsx'],[EOdataloc,[eegdir, '.out']]);  % takes an excel file (of votes) as input and creates a similar .out file for reconstructing the EEG data (P)
        % Just see how the file is, I created 'CPT00112_forSOBI_voteforrecon.out' by using 'createOutFromXLS' (P)
        
        outfile = [eegdir, '.out']; % 'CPT03811_forSOBI.out' This is the file created by createOutFromXLS (P)
%         CreateNoiseReconFile([SmartOutloc, nameoffile, '_voteforrecon_EO.xlsx'],[EOdataloc,[eegdir, '_EyesOpenNoise.out']]);
%         noiseoutfile = [eegdir, '_EyesOpenNoise.out'];
%         %write the reconstructed data
       reconEEGdata_forCPTdata(eegdir, 1, EOdataloc,  outfile); % ['CPT03811_forSOBI',1,'\\Dss02544-cmb-d\h\CPT_SOBI_Output\CPT03811_forSOBI\CPT03811_forSOBI_Sources_CPT\','CPT03811_forSOBI.out'  ]  (P)
       % writes clean reconstructed data in terms of .dat, .generic and .evt files
       
       reconEEGdata_forCPTdata(eegdir, 2, EOdataloc,  outfile); % ['CPT03811_forSOBI',2,'\\Dss02544-cmb-d\h\CPT_SOBI_Output\CPT03811_forSOBI\CPT03811_forSOBI_Sources_CPT\','CPT03811_forSOBI.out'  ]  (P)
       % writes noise data (also called reconInverse) in terms of .dat, .generic and .evt files
       
       
       %  reconNoise_forRestdata(eegdir, 1, EOdataloc, outfile, noiseoutfile);
       % reconNoise_forRestdata(eegdir, 2, EOdataloc, outfile, noiseoutfile);
        IntegrateResponseTrigs(eegfile);
    end
end


function IntegrateResponseTrigs(nameoffile)

    reconevtfile = strcat('recon_',nameoffile, '.evt'); %original              This is the event file which got created in the file reconEEGdata_forCPTdata.m (P)
    %reconevtfile = strcat('recon_',nameoffile, '_CPT.evt'); %changed for the test data. change back to orig for normal files
    %eegfile = strcat(nameoffile, '_epoched_forSOBI'); %original
    eegfile = nameoffile; %changed for the test data. change back to orig for normal files
    resptrigevtfile = strcat(eegfile(1:8), '_GoodTrials.evt');    % removed 'for_SOBI' and added '_GoodTrials.evt' (P)
    outloc = '\\Dss02544-cmb-d\h\CPT_SOBI_Output\';   
    recondataloc = strcat(outloc, nameoffile,'\', [nameoffile, '_Sources_CPT\'], [nameoffile, '_recon\']); % folder for reconstructed data (P)
    Trigdataloc = '\\Dss02544-cmb-d\h\CPT_SOBI_data\';
        
    %read in the timing and the triggers of the reconstructed data output
    reconTimes = textread(strcat(recondataloc,reconevtfile),'%s');   % everything falls as one single column. So basically, to access a certain parameter like 'tmu', we need to jump by the index of 4 because there were 4 colums in the event file (P)
    recontmu = [];
    recontrig = [];
    for i =5:4:size(reconTimes) % 5 is the starting point because 1st row (4 elements) was header (P)
        recontmu(end+1) = str2num(cell2mat(reconTimes(i))); % copying time stamps (P)
    end

    for i =7:4:size(reconTimes)  % 7 is the starting point because 1st row (4 elements) was header (P)
        recontrig(end+1) = str2num(cell2mat(reconTimes(i))); % copying triggers (P)
    end

    %Read in the file with the response triggers included.
    RespTrigTimes = textread(strcat(Trigdataloc,resptrigevtfile), '%s'); % This is the preSOBI GoodTrial event file which has the information about the original response triggers timings (P)
    % everything falls as one single column. So basically, to access a certain parameter like 'tmu', we need to jump by the index of 9 because there were 9 colums in the event file (P)
    
    faRTs = [];  %hahaha, "farts"... false alarm response times
    HitRTs = []; %hit response times 
    nontarget = 0; %long line / trigger 2
    target = 0; %short line / trigger 4
    hit = 0; %response after a short line / trigger 4 then 1
    falsealarm = 0; %response after a long line/ trigger 2 then 1
    corrrej = 0; %no response after a long line / trigger 2 then 2 or 4
    miss = 0; %short line and no response / trigger 4 then trigger 2 or 4
    faPosition = []; %position of each false alarm in the event series
    hitPosition = []; % position of each hit in the event series
    NumStims = 0; %index current position

    for i = 13:9:size(RespTrigTimes)-4  % I don't think -4 would have any effect.(P)
    
        if strcmp(RespTrigTimes(i),'2')
            NumStims = NumStims+1;    
            nontarget = nontarget + 1;
        elseif strcmp(RespTrigTimes(i),'4')
            NumStims = NumStims+1;
            target = target + 1;  
        end
    
        if strcmp(RespTrigTimes(i),'1') && strcmp(RespTrigTimes(i-9),'2')
            falsealarm = falsealarm +1;
            faResponseRT = cell2mat(RespTrigTimes(i-4)); %record the absolute RT            % This is the time stamp of response trigger (P)
            nonTargetTime = cell2mat(RespTrigTimes(i-13));%record the target presentation period        % This is the time response of the corresponding stimulus trigger (P)
            faRTs(end+1) =  str2double(faResponseRT) - str2double(nonTargetTime); %calculate the relative RT         % They could have simply used the RT mentioned in event file (P)
            faPosition(end+1) = NumStims; %record the position of the FA within the event series 
        
        elseif strcmp(RespTrigTimes(i),'1') && strcmp(RespTrigTimes(i-9),'4')
            ResponseRT = cell2mat(RespTrigTimes(i-4));
            TargetTime = cell2mat(RespTrigTimes(i-13));
            HitRTs(end+1) =  str2double(ResponseRT) - str2double(TargetTime);
            hitPosition(end+1) = NumStims; %record the position of the hit within the event series
            hit = hit + 1;
        
        end

        if strcmp(RespTrigTimes(i),'2') && strcmp(RespTrigTimes(i-9),'2')
            corrrej = corrrej +1; %if there's no response after the '2' trigger, it's a correct rejection.
        elseif strcmp(RespTrigTimes(i),'2') && strcmp(RespTrigTimes(i-9),'4')           
                miss = miss + 1; %if there's no response after the '4' trigger, it's a miss.
        end

    end

    writeBESAevt2(strcat(recondataloc,'recon_',eegfile,'_withRespTrigs.evt'), recontmu,recontrig,HitRTs, faRTs,faPosition, hitPosition, NumStims);  % Integrates response triggers but not with the actual timestamp (P)
    disp('Response Triggers successfully integrated.');
end              

function writeBESAevt2(evtfile, recontmu, recontrig, HitRTs,faRTs, faPosition, hitPosition, NumStims)
    % creating event file specifically for CPT epoched data.
    fid = fopen(evtfile,'w');
    if fid == -1
        fprintf(1,'Error creating event file with integrated response triggers.\n');
        return;
    end
    
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    targethitcounter = 1;
    nontargetfacounter = 1;
    
    for t = 1:1:size(recontrig,2)  % these are the number of rows in the event file formed during reconstruction
        fprintf(fid, '%d\t%d\t%d\t%s\n', recontmu(t)+1.5e6, 1, recontrig(t), strcat('Trial:',num2str(t)));   % 1.5 s added to time stamps because the epoch created had 1.5 sec as pre-stim interval
        if strcmp(num2str(recontrig(t)),'2')
            if faPosition ~=0
                if strcmp(num2str(t),num2str(faPosition(nontargetfacounter)))
                    fprintf(fid, '%d\t%d\t%d\t%s\n', recontmu(t)+2e6, 1, 1, strcat('Trial:',num2str(t)));   % Why was 2 sec added arbitrarily to the time stamp? (P)
                    if nontargetfacounter < size(faPosition,2)
                        nontargetfacounter = nontargetfacounter + 1;
                    end
              
                end  
            end
        end
           
        if strcmp(num2str(recontrig(t)) ,'4')
            if hitPosition ~=0
                if strcmp(num2str(t),num2str(hitPosition(targethitcounter)))
                    fprintf(fid, '%d\t%d\t%d\t%s\n', recontmu(t)+2e6, 1, 1, strcat('Trial:',num2str(t))); %time added to recontmu is based on prestimulus period + set time for response trigger just for marking the epoch with a response trigger -- it is not set at the time of response.
               
                    if targethitcounter < size(hitPosition,2)
                        targethitcounter = targethitcounter + 1;
                    end
                end
            end
        end
    end
    fclose(fid);  
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