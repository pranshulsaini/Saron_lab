% This code is developed to run SOBI on EEG data; 
% author: Manish Saggar (mishu@cs.utexas.edu). 
% Iman has changed the code and added the embedded source trial  ver:04/12/2012
% (irezazadeh@ucdavis.edu)
%Chivon has changed the code to run with STR task data. 

% Last modified before Pranshul: 18th Feb 2015


function y = tSOBI_STR_betterSMART(eegfilename, eegfile, evtfile, sfpfile, outputDir, chunks)
 
    diary tSOBI_logFile;   % writes everything printed on the command window
    eegfile = strtrim(eegfile);  % removes leading and trailing whitespace characters from txt       
    
    % Variables
    hdr = readBESAsb_header(eegfile); % this returns a structure having all the information which a generic file has (P)
    tang_tau = [1:10,12:2:20,25:5:100,120:20:300]; % in msec                These are the times for which we will check correlation (P)

    
    % derive tau by changing tang_tau according to our Fs
    tau = [0 unique(round(hdr.Fs*(tang_tau./1000)))]; % this is like the number of sample points for the time points mentioned in tang_tau (P)

    sinThreshold = 1e-08;         % Threshold for non-diagonal extent during SOBI (P)
 
    
   % Creating trial of interest (toi) variable by reading data from the Event file.
    evts = textread(evtfile,'%s');
    tmu = []; 
    trig = [];
     
%     for i = 5:1:size(evts,1)
%         if mod(i,5) == 0
%             tmu(end+1) = str2num(cell2mat(evts(i)));
%         elseif mod(i,5) == 2
%             trig(end+1) = str2num(cell2mat(evts(i)));
%         end
%     end 
    %This is code to read in the goodtrials file so the event file doesn't
    %have to be manually created for every file.
    % %original info
    origtrig = [];
    origtmu = [];
    %origISI = [];    % I don't have this information in stroop event file. So commented (P)
    %origBlock = [];    % I don't have this information in stroop event file. So commented (P)
    %origRT = [];  % I don't have this information in stroop event file. So commented (P)
    %origCount = [];  % I don't have this information in stroop event file. So commented (P)

    for i = 5:7:size(evts,1)    % the gap of 7 is there because evts file has everything in one column. So, all the rows are converted to a column (P)
        origtmu(end+1) = str2num(cell2mat(evts(i)));    % at the end of the loop, the array are of dimension: 1x1032 (P)
        origtrig(end+1)= str2num(cell2mat(evts(i+2)));  % at the end of the loop, the array are of dimension: 1x1032 (P)
        %origISI(end+1) = str2num(cell2mat(evts(i+6)));   % at the end of the loop, the array are of dimension: 1x1032 (P)
        %origBlock(end+1) = str2num(cell2mat(evts(i+5)));   % at the end of the loop, the array are of dimension: 1x1032 (P)
        %origRT(end+1) = str2num(cell2mat(evts(i+7)));   % at the end of the loop, the array are of dimension: 1x1032 (P)
        %origCount(end + 1) = str2num(cell2mat(evts(i+8)));   % at the end of the loop, the array are of dimension: 1x1032 (P)
    end

    Trialstep = [];
    %I want to drop the response triggers to make SOBI process a little less chaotic
    
    for i = 1:1:size(origtrig,2)    % size(origtrig,2)  = 1032 (P)
        if origtrig(i) >1000 % This will let pass only the stimulus triggers
            trig(end+1) = origtrig(i);    % so this is the trigger for stimuli (P)
            tmu(end+1) = origtmu(i);  % just the corresponding timestamp (P)
        end
    end
     

    toiAll = [];
    trigAll = trig;
    cAll = 1;
    tmu = tmu./1e6; % converting to seconds to multiply by # of samples.
    for i = 1:1:size(trig,2)
        toiAll(cAll,1) = round(tmu(i)*2048) - 512; %Set the stim-locked epoch range in samples            512 ~ 250 ms      (P)
        toiAll(cAll,2) = round(tmu(i)*2048) + 3584;  % 3584 ~ 1750 ms (P)
        cAll = cAll + 1;
    end
    
      
    % dividing the toiAll into chunks based upon the value of parameter
    % chunks. This parameter is always set to 1 for the STR data.
    
  
    if chunks == 1  % As chunk is always 1, only this condition will be implemented

        % calc correlation matrices
        Rx = createCorrMat(hdr, eegfile, toiAll, tau); % returns a matrix of 88 x 88 x n_tau

        % Feed Correlation Matrices to Joint Diagonalization.
        [W, B, D] = jointDiagonalization(Rx, sinThreshold);  % If we multiply with W with the input data, we will get independent components. D (m x n*m)is the collection of almost diagonl matrices (P). B is pre-sphering matrix

        W_scaled = real(scaleW(W));     % I do not know why it is required (P)
        size(W_scaled)
        % output nSources x nTrials per condition; S = WX;            
        aST = createSourceTrials(hdr, W_scaled, eegfile, toiAll, 'STR', eegfilename, trig); %m x n_trial*epoch_time

        % averaged source components across all epoch;             
        eST = createEvokedSourceTrial(hdr, W_scaled, eegfile, toiAll, 'STR', eegfilename);  % m x 4096 (P)

        % averages source components across condition/trigger specific epochs  (P)
        embeddedST = createEmbeddedEvokedSourceTrial(hdr, W_scaled, eegfile, toiAll, 'STR',evtfile,eegfilename,trig,tmu); % m x (4096*5)

        % create HTML file using SMART
        smart_autism(aST, eST, hdr, W_scaled, eegfilename, toiAll, trigAll, 'STR', sfpfile);
        % The above command writes aST, eST, new EOG and EMG artifact information, new peaks around power noise,... all that in the HTML file (P)
        
    else               % I have not gone through this because I think we will also keep chunks =1 for Stroop data
        totTrials = size(toiAll,1);
        toiAll_trls = {};
        trigAll_trgs = {};
        toiAll_ind = 1;
        for i = 1:1:chunks
            if i ~= chunks
                toiAll_trls{i} = toiAll([toiAll_ind:toiAll_ind+floor(totTrials/chunks)-1],:);
                trigAll_trgs{i} = trigAll(1,[toiAll_ind:toiAll_ind+floor(totTrials/chunks)-1]);
                toiAll_ind = toiAll_ind + floor(totTrials/chunks);
            else
                toiAll_trls{i} = toiAll([toiAll_ind:end],:);
                trigAll_trgs{i} = trigAll(1,[toiAll_ind:end]);
                toiAll_ind = toiAll_ind + floor(totTrials/chunks);
            end
        end
        
        for i = 1:1:chunks
            % running SOBI once on all conditions together
            % calc correlation matrices
            Rx = createCorrMat(hdr, eegfile, toiAll_trls{i}, tau);

            % Feed Correlation Matrices to Joint Diagonalization.
            [W, B, D] = jointDiagonalization(Rx, sinThreshold);
            W_scaled = real(scaleW(W));    

            % output nSources x nTrials per condition; S = WX;            
            aST = createSourceTrials(hdr, W_scaled, eegfile, toiAll_trls{i}, strcat('ACC',num2str(i)));

            % output projection of each source in sensor space for the 
            % averaged epoch;            
           % eST = createEvokedSourceTrial(hdr, W_scaled, eegfile, toiAll_trls{i}, strcat('ACC',num2str(i)));
            %eSTembedded = createEmbeddedEvokedSourceTrial(hdr, W_scaled, eegfile, toiAll_trls{i}, strcat('ACC',num2str(i)));

            % create HTML file using SMART
            smart_autism(aST, eST, hdr, W_scaled, eegfile, toiAll_trls{i}, trigAll_trgs{i}, strcat('ACC',num2str(i)), sfpfile);
        end
    end
    diary off;             
end


function smart_autism(aST, eST, hdr, W, eegfile, toi, trig, cond, sfpfile)  % [aST, eST, hdr, W_scaled, eegfilename, toiAll, trigAll, 'STR', sfpfile]
    emgThreshold = 0.6;   % How are these thresholds decided? Is this what Cliff was talking about as 'classifier by Manish'?
    eogThreshold = 0.7;   % How are these thresholds decided? Is this what Cliff was talking about as 'classifier by Manish'?
   
    update_sfp(sfpfile);   % removes the unnecessary rows in the beginning
    
    [tmp1, tmp2, tmp3] = fileparts(sfpfile); 
    newsfpfile = strcat(tmp1,'\',tmp2,'_updated',tmp3);
   
    el = readlocs(newsfpfile);  % reads electrode location coordinates and other information from a file. (P)
    %output is a structure containing the channel names and locations (if present). It has three fields: 'eloc.labels', 'eloc.theta' and 'eloc.radius' (P)
    % I get this warning: coordinate conversion failed. However, the output seems fine. It has 5 columns: [electrode number, labels, X, Y, Z] (P)
    
    % calculate EMG components based on Auto-correlation function on the evoked (average) source trial.                ACF means autocorrelation function (P)

    [emgnew, ac] = findACFbasedEMG(aST, W, toi, emgThreshold);   % ac (dim: m X 64) contains mean autocorrelation coeff, and emgnew contains the list of found emg channels (P)
    
    
    % calculate EOG components based on Auto-correlation function on the evoked (average) source trial.          Electrooculography (EOG) is a
    % technique for measuring the corneo-retinal standing potential that exists between the front and the back of the human eye (P)
    eognew = findACFbasedEOG(ac, W, eogThreshold);  % eog new contains the list of found EOG sources (P)
    
    % calculate Peak components based on PowerSpectrum analysis on each Source trial.                  This is kind of finding the bad peaks for power noise around 60 Hz (P)
    [peaksnew peaks Sps freq]= findPeaks(aST, hdr, W, toi);
    % peaksnew will have the sources indexes of the peaks (P)
    % peaks will have the values of the frequencies at which peak happens (P)
    % Sps will have the mean power spectrum of all sources. Dim: 88 x nFFT (P)
    % freq will have the list of frequencies for which we have the spectra (P)
    
    [tmp1, tmp2, tmp3] = fileparts(eegfile); 
    eegfile = tmp2; % without full path
    printHTML(emgnew, eognew, peaksnew, peaks, ac, el, W, Sps, aST, eST,...
        eegfile, cond, freq, toi); 
    
    save(strcat(strrep(eegfile,'.dat',''),'_',cond,'_W'),'W','aST','toi','trig');   % It will create the file and store 'W','aST','toi','trig' in that (P)
end


function [y acSt] = findACFbasedEMG(St, W, toi, thresh)  % ACF means autocorrelation function (P)
    y = [];
    acSt = [];
    blockSize = 4096;
    for s = 1:1:size(W,1)  % size(W,1) = m
        acst = []; 
        
        for t = 1:1:size(toi,1)  % n_trials = 256
            %  ac = acf_old(St(s, toi(t,1):toi(t,2))',64);
             ac = acf_old(St(s, (t-1)*blockSize+1:t*blockSize)',64);  % Estimate the coefficients of the autocorrelation (covariance with its own lagged value). 64 is the max lag considered in digital values (P)
             acst(t,:) = ac.ac; % autocorrelation coefficients. Final size of acst would be: 256 x 64   (P)
        end

        % average the acf for all trials for each source.
        acSt(s,:) = mean(acst,1);    % Size of acSt for all sources: 88 x 64

        if min(acSt(s,1:20)) < thresh  % trying min instead of mean             Manish Saggar's paper figure shows that EMG autocorrelation decreases sharply with lag, unlike neural signals (P)
            y = [y s];          % source index is added as an emg channel (P)
        end
    end
end


function [y] = findACFbasedEOG(ac, W, thresh)
    y = [];
    for i = 1:1:size(W,1)
       if mean(ac(i,:)) > thresh   % here mean value (P)
           y = [y i];
       end
    end
end


function [y pp psStNL freq] = findPeaks(St, hdr, W, toi)
    y =[];
    pp = [];
    NFFT = 2700;    % controls the frequency resolution.             Thi is being used as segment length, not NFFT (P)        
    % NFFT specifies the number of discrete Fourier transform (DFT) points to use in the PSD estimate. The default nfft is the greater of 256 or the
    % next power of 2 greater than the length of the segments (P). 
    
    %Hs = spectrum.welch({'Hann','periodic'}, NFFT, 0);  % NFFT act as segment length here. 0 overlap (P)
    psSt = [];
    psStNL = [];
    blockSize = 4096;
    
    % calculate power spectrum for each trial and each source.
    for s = 1:1:size(W,1)   % 88 (P)
        psst = [];
        for t = 1:1:size(toi,1)  % n_trials = 960
           % tmp = St(s, toi(t,1):toi(t,2));
            tmp = St(s, (t-1)*blockSize+1:t*blockSize);
            %ps = psd(Hs, tmp,'Fs',hdr.Fs);          % I will have to use pwelch() because psd() is deprecated. ps should be of dim: (2700/2 + 1) x 1 = 1351 x 1
            
            % I will have to use one of the following commands. Both gives the same answer. NFFT is being used as window length as well as the NFFT (P)
            %[pxx, f] = pwelch(a,hann(NFFT,'periodic'),0, NFFT,2048 ); (P)
            [pxx, f] = pwelch(tmp,hanning(NFFT,'periodic'),0, NFFT,2048 ); %(P)
            %The length of ps would be  (nFFT/2 + 1) (P)
            
            psst(t,:) = pxx';  % I would have to write psst(t,:) = pxx'. In the end, psst will have dim: 256 x 1351; (P)
        end
        
        % average the power spectrum for all trials for each source.
        psSt(s,:) = 10*log10(mean(psst,1));   % In the end, psSt will have the dimension: m x 1351 (P)
        psStNL(s,:) = mean(psst,1); % In the end, psSt will have the dimension: m x 1351 (P)
        freq = f; % I would have to write freq = f; (P)

        % to avoid error from findpeaks                    
        if max(psStNL(s,:)) <.5,                 
            continue;                   % It means it will skip all the commands below and enter next for-loop for next value of s
        end

        [p,l] = findpeaks(psStNL(s,:),'SORTSTR','descend','MINPEAKHEIGHT',.5); % sortstr = descend means that the largest peak will get number 1, the second largest 2, and so on (P)
        % p will carry the peak values. l will carry the corresponding indexes
        
        npeaks = 5;  % these are the maximum peaks we will allow (P)
        
        if length(l) < 1        % no peaks identified, which means everything was less than 0.5  (P)
            continue;
        elseif length(l) < 5    
            npeaks = length(l);
        end

        freqPeaks = round(freq(l));
        
        for i = 1:1:npeaks   % just look at first three peaks            It should be "five" peaks in the comment (P)
           if freqPeaks(i) == 59 || freqPeaks(i) == 60 || freqPeaks(i) == 61          % I think this is for power noise (P)
                y = [y s];  % The source index is added. In the end, it length would be (n_59or60or61_peaks x 1) (P)
                pp = [pp freqPeaks(i)];  % Peaks are added. In the end, it length would be (n_59or60or61_peaks x 1) (P)
           end
        end
        
    end
end


function sT = createEvokedSourceTrial(hdr, W, eegfile, toi, cond, eegfilename) % hdr, W_scaled, eegfile, toiAll, 'STR', eegfilename (P)
    blockSize = 4096;
    sT = zeros(size(W,1), blockSize);   % 88 x 4096 (P)
    for t = 1:1:size(toi,1)
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2));
        tmp = W*x;
        % demeaning the sources from -80: 20  msec (this is shifted by 1 sample also).
        tmp = tmp - mean(tmp(:,348:553),2)*ones(1,size(tmp,2));
        sT = sT + tmp;   % sT is summing up source projection fo all blocks/epochs
    end
    
    sT = sT./size(toi,1);  % scaling down because it contains the values added across all the trials (P)
%     sT = sT * 1e3;
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path    
    writeBESAsb_data(sT, strcat(eegfilename,'_',cond,'_EvokedSourceTrial.dat'));   % writing averaged evoked data (P)
    
    % creating event file
    fid = fopen(strcat(eegfilename,'_',cond,'_EvokedSourceTrial.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Sensor Projection\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
  %  for t = 1:1:size(W,1)
       fprintf(fid, '%d\t%d\t%d\t%s\n', 250*1e3, 1, 1, 'Stimulus');   % 250 is chosen because epoched data has 250 ms pre stimulus time (P)
  %  end
    fclose(fid); 
    
    % creating evoked source trial projections
    chData = zeros(size(W,1), size(W,1)*blockSize);   % m x (m*4096) (P)
    for t = 1:1:size(W,1)  % for each source
        A = inv(W);
        chData(:,(t-1)*blockSize+1:t*blockSize) = A(:,t)*sT(t,:); % append projection of each source into sensor space
    end
    writeBESAsb_data(chData, strcat(eegfilename,'_',cond,'_EvokedSourceTrialSensorProj.dat'));
    
 % creating event file
    fid = fopen(strcat(strrep(eegfile,'.dat',''),'_',cond,'_EvokedSourceTrialSensorProj.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Sensor Projection\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    for t = 1:1:size(W,1)  % m times. The loop is to merely say that this is the same for all sources. See that the last argument has t (P)
       fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048))*1e6+1, 1, t, strcat('Trig.',num2str(t))); % Every row corresponds to projection of the corresponding source
    end
    fclose(fid);   
       
end


function sT = createEmbeddedEvokedSourceTrial(hdr, W, eegfile, toi, cond,evtfile, eegfilename,trig,tmu) % [hdr, W_scaled, eegfile, toiAll, 'STR',evtfile,eegfilename,trig,tmu]
    blockSize = 4096;
    no_event=5;  % these are number of different event conditions, i.e. {incongruent, incongruent_neutral, neutral, congruent_neutral, congruent}
    sT_new = zeros(size(W,1), blockSize,no_event); % mx 4096 x 2 (P)
    %sT = zeros(size(W,1), blockSize*no_event);
    sT = [];
    trigger_counter=zeros(no_event,1); % 2x2. It was zeros(no_event). It should be zeros(no_event,1). I fixed it (P)
    
%     evts1 = textread(evtfile,'%s');
%     tmu = [];
%     trig = [];
%     for i = 5:1:size(evts1,1)
%         if mod(i,5) == 0
%             tmu(end+1) = str2num(cell2mat(evts1(i)));
%         elseif mod(i,5) == 2
%             trig(end+1) = str2num(cell2mat(evts1(i)));
%         end
%     end 
%     
    for t = 1:1:size(toi,1)
        trigger_no=trig(t);  % it would be a number
        trigger_no = int2str(trigger_no);
        condition = str2double(trigger_no(1)); % Ut would be either 1, 2, 3, 4, 5
        
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2));
        tmp = W*x;
        
        % demeaning the sources from -80: 20  msec (this is shifted by 1 sample also).
        tmp = tmp - mean(tmp(:,348:553),2)*ones(1,size(tmp,2));
        sT_new(:,:,condition) = sT_new(:,:,condition) + tmp; % adding them separately for the 5 condition
        trigger_counter(condition)=trigger_counter(condition)+1;  
    end
    
    %sT = sT./size(toi,1);
    for i=1:no_event  % normalisation (P)
        sT_new(:,:,i) = sT_new(:,:,i)./trigger_counter(i);
        sT=[sT sT_new(:,:,i)];       % In the end, it would be m x 4096 x 5
    end
    
%     sT = sT * 1e3;
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path    
    writeBESAsb_data(sT, strcat(eegfilename,'_',cond,'_EmbeddedEvokedSourceTrial.dat'));  % wrote the conditon/trigger specific evoked source potential (P)
    
    % creating evoked source trial projections
    
%    pData = zeros(size(W,1), size(W,1)*blockSize,event_no);
    chData = zeros(size(W,1), size(W,1)*blockSize*no_event);    % 88x (88*4096*5)
    bias_time=0;   
    for t = 1:1:size(W,1)  % for each source
        for i=0:no_event-1
            % append projection of each source into sensor space
            A = inv(W);  % could be put outside the loop (P)
            %chData(:,(t+i-1)*blockSize+1+bias_time:(t+i)*blockSize+bias_time) = A(:,t)*sT_new(t,:,i+1);
            chData(:,(i)*blockSize+1+bias_time:(i+1)*blockSize+bias_time) = A(:,t)*sT_new(t,:,i+1);  % reconstructing the signal (P)
        end
        bias_time=((1+i)*blockSize+bias_time);
    end 
   
    writeBESAsb_data(chData, strcat(eegfilename,'_',cond,'_EmbeddedEvokedSourceTrialSensorProj.dat')); % wrote the conditon/trigger specific evoked sensor potential (P)
        
    % creating event files: one with conditions and another with sources
    fid = fopen(strcat(eegfilename,'_',cond,'_EmbeddedEvokedSourceTrialSensorProj.evt'),'w');
    fid2 = fopen(strcat(eegfilename,'_',cond,'_EmbeddedEvokedSourceTrialSensorProjSourceNums.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Embedded Sensor Projection\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    fprintf(fid2,'Tmu\tCode\tTriNo\tComnt\n');
    bias_time=0;
    for t = 1:1:size(W,1)
       for i=0:no_event-1
         %   fprintf(fid2, '%d\t%d\t%d\t%s\n', ((t+i-1)*blockSize+1)*1e3+bias_time, 1, t, strcat('Trig.',num2str(t)));
             fprintf(fid2, '%d\t%d\t%d\t%s\n', ((i)*(blockSize/2048)*1e6+bias_time)+1, 1, t, strcat('Trig.',num2str(t))); % with sources (P)
            fprintf(fid, '%d\t%d\t%d\t%s\n', ((i)*(blockSize/2048)*1e6+bias_time)+1, 1, (i+1), strcat('Trig.',num2str(t))); % with conditions (P)
       end 
       bias_time=((i+1)*(blockSize/2048))*1e6+bias_time;
    end
%     for t = 1:1:size(trig,2)
%        fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048))*1e6+1, 1, trig(t), strcat('Trial:',num2str(t)));
%     end
    fclose(fid);
    fclose(fid2);    
end


function sT = createSourceTrials(hdr, W, eegfile, toi, cond, eegfilename, trig)   % [ hdr, W_scaled, eegfile, toiAll, 'STR', eegfilename, trig] (P)
    subjid = eegfilename;    % eegfilename does not have 'dat' in it (P)
    output = strcat(subjid,'_Sources_',cond);  % eg:  STR00412_Sources_STR
    
    if exist(output,'dir') ~= 7
       fprintf(1, '\nArtifacts dir doesn''t exist, Creating one...');
       mkdir(output);
    else 
       fprintf(1, '\nArtifacts dir already exists, Overwriting...'); 
       mkdir(output);
    end
    
    cd(output); % changes the current path of matlab to 'output' directory    
    blockSize = toi(1,2)- toi(1,1);
    sT = zeros(size(W,1), size(toi,1)*(toi(1,2)- toi(1,1))); % m x n_trial*epoch_time
    
    for t = 1:1:size(toi,1)
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2)); % 88 x n_points_for_1_epoch (P)
        tmp = W*x;    % independent source signal. 88 x n_points_for_1_epoch. So, SOBI is done separately for each epoch (P)
        
        % demeaning the sources from -80: 20 msec (shifted by 1 sample).
        tmp = tmp - mean(tmp(:,348:553),2)*ones(1,size(tmp,2)); % As the my calculation, the time points for 1 trial are 4096, according;y [348, 553].
        sT(:, (t-1)*blockSize+1:t*blockSize) = tmp;  % storing independent source signal for all epochs (P) 
    end

% I am commenting it (Pranshul). I think it wastes a lot of memory and processing time as well
%     A = inv(W);  % This is required for reconstruction (P)
%     [tmp1, tmp2, tmp3] = fileparts(eegfile); % [filepath,name,ext] = fileparts(file) (P)
%      
%     eegfile = tmp2; % without full path  
%     clear tmp tmp1 tmp2 tmp3;
%     for s = 1:1:size(W,1)  % for each source
%         chData_tmp = zeros(size(W,1), size(toi,1)*blockSize);    %88 x n_trial*epoch_time (P)
%         for t = 1:1:size(toi,1)
%               
%            chData_tmp(:,(t-1)*blockSize+1:t*blockSize)= A(:,s)*sT(s,(t-1)*blockSize+1:t*blockSize);  % reconstructing original signal (P)
%            %chData_tmp(:,(t-1)*blockSize+1:(t*blockSize)-(blockSize/2))= A(:,s)*sT(s,(t-1)*blockSize+1:(t*blockSize)-(blockSize/2));
%            %chData_tmp(:,(t*blockSize)-(blockSize/2)+1:t*blockSize)= A(:,s)*sT(s,(t*blockSize)-(blockSize/2)+1:t*blockSize);
%            
%         end 
%         writeBESAsb_data(chData_tmp, strcat(eegfilename,'_',cond,'_Source_',num2str(s),'_Projection_Into_SensorSpace_AllTrials.dat'));    % writing reconstructed data file for each source (P)
%         clear chData_tmp   % I don't think it is required since it has been assigned zero in the beginning of the loop (P)
%     end 
%     
%     % creating event file for sensorspace projection files
%     fid = fopen(strcat(eegfilename,'_',cond,'_SensorSpace_AllTrials.evt'),'w');  % timestamp will not be actual timestamp
%     fid2 = fopen(strcat(eegfilename,'_',cond,'_SensorSpace_AllTrialsTrigs.evt'),'w'); % timestamp will not be actual timestamp
%     
%     if fid == -1
%         fprintf(1,'Error creating Event File for Source Trials\n');
%     end
%     
%     fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
%     fprintf(fid2,'Tmu\tCode\tTriNo\tComnt\n');
%     
%     for t = 1:1:size(toi,1)   % I didn't understand the need to create these files (P)
%          %fprintf(fid, '%d\t%d\t%d\t%s\n', toi(t,1)*1e3, 1, t, strcat('Trial:',num2str(t)));
%         fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048)*1e6+1), 1, t, strcat('Trial:',num2str(t)));    % this line has t. The tmu here will not be actual tmu   (P)
%         fprintf(fid2, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048)*1e6+1), 1, trig(t), strcat('Trial:',num2str(t)));  % this line has trig(t). The tmu here will not be actual tmu   (P)
%     end
%     
%     fclose(fid);
%     fclose(fid2);

    % writing the source space data into files
    %     sT = sT * 1e3;    
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path    
    writeBESAsb_data(sT, strcat(eegfilename,'_',cond,'_SourceTrials.dat'));   %writes .dat and .generic files for the source_space data

    % creating event file
    fid = fopen(strcat(eegfilename,'_',cond,'_SourceTrials.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Source Trials\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    for t = 1:1:size(toi,1)
         %fprintf(fid, '%d\t%d\t%d\t%s\n', toi(t,1)*1e3, 1, t, strcat('Trial:',num2str(t)));
        fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048)*1e6+1), 1, t, strcat('Trial:',num2str(t))); % this is the same one as sensor projection event file (P)
    end
    fclose(fid);
end


function writeBESAsb_data(X_cap, new_file)  % [sT, strcat(eegfilename,'_',cond,'_SourceTrials.dat')] (P)
     fid = fopen(new_file, 'w','ieee-le');
    if fid == -1
       fprintf(1,'In writeBESAsb_data:: Output file not found\n');
       exit;
    end
    fwrite(fid, X_cap, 'float32');  
    [hdr.nChans,hdr.nSamples] = size(X_cap);
    hdr.Fs = 2048;
    writeBESAgeneric(new_file,hdr);
    fclose(fid);
end


function writeBESAgeneric(new_file, hdr)
%     [path, name, ext, ver] = fileparts(new_file);
    fid = fopen(strrep(new_file,'.dat','.generic'), 'w');
    fprintf(fid,'BESA Generic Data\n\n');
    fprintf(fid,'nChannels=%i\n\n', hdr.nChans);
    fprintf(fid,'sRate=%f\n\n', hdr.Fs);
    fprintf(fid,'nSamples=%i\n\n', hdr.nSamples);
    fprintf(fid,'format=float\n\n');
    fprintf(fid,'file=%s\n\n', new_file);
    fclose(fid);
end


% tau values should be in samples. 
function Rx = createCorrMat(hdr, eegfile, toi, tau)   % toi contains the time ranges [pre_stim post_stim]
    Rx = zeros(hdr.nChans, hdr.nChans, length(tau));  % tau contains the index time points for which correlation has to be calculated (P)
    %blockSize = 2999;
    % Read data. We are dividing continuous data into epochs in order to use only useful data to obtain matrices to perform SOBI (P)
    for t = 1:1:size(toi,1)    % size(toi,1) = 1031. The total number of trials are 960. The surplus is due to the addition of one-rows in good trials
        % the data should be read in nChans x Time format
        % x = readBESAsb_data(eegfile, hdr, toi(t,1)*2+1, toi(t,2)*2 + max(tau));
        x = readBESAsb_data(eegfile, hdr, toi(t,1), toi(t,2) + max(tau));    % Returns a 88 x n_time_points matrix. max tau is added to incorporate the max temporal correlation we want (P)
        
        % de-mean data
        if isnan(x)~=1
        x = x - mean(x,2)*ones(1, size(x,2));
        
        % calculating corr matrix at different tau values
        for i = 1:1:length(tau)
          %  tt=size(x); 
          %  if tt(1,2)==toi 
                
               Rxx = x(:,1:toi(t,2)-toi(t,1)) * x(:,1+tau(i):toi(t,2)- toi(t,1)+tau(i))';  % 88 x 88 matrix. They don't divide by the variance because it is not mentioned in SOBI paper (P)
              % Rxx = x(:,1:blockSize) * x(:,1+tau(i):blockSize+tau(i))'; 
               Rx(:,:,i) = Rx(:,:,i) + 0.5*(Rxx + Rxx');   % I think this is done because of the SOBI formulation. However, Rx(:,:,i) on RHS shows that correlation b/w every channel for each tau is added every epoch
%                if isnan(Rxx)==1
%                        disp(strcat (i,'>>>>',t));
%                        pause;
%                end
%                
%                 if isnan(Rx)==1
%                        disp(strcat (i,'>>>>',t));
%                        pause;
%                end

           % end 
        end
        fprintf(1,'\n Calculating correlation matrices, for epoch: %d', t);
        end
     end
    Rx = Rx./size(toi,1);  % This is done because we Rx is added for every epoch as it is looped for each epoch, i.e. size(toi,1)
end


function [W, B, D] = jointDiagonalization(Rx, sinThreshold)
    
    % calculating pre-sphering matrix
    B = preWhitenData(Rx(:,:,1)); % m x m, where m = number of channels in the EEG data (P)
    
    % whitening the correlation matrices
    Rx_white = zeros(size(Rx,1),size(Rx,2),size(Rx,3)-1);  % Dim: m x m x (n_tau-1). leaving the instantaneuous correlation.   Why??? (P)
    
    for i = 1:1:size(Rx,3)-1
      Rx_white(:,:,i) = B * Rx(:,:,i+1) * B';   % Why are not using the corresponding B, but the fixed one?? (P)
    end    
    
    % calculating joint diagolization matrix
    fprintf(1,'\n====Joint Diagonalization of Correlation Matrices====\n');
    [W, D] = joint_diag(reshape(Rx_white,size(Rx_white,1),size(Rx_white,2)*size(Rx_white,3)), sinThreshold); % W is an m*m unitary matrix. D (m x n*m)is the collection of almost diagonl matrices (P) 
    W = W' * B;              % if we multiply this W (m x m) with the original signal, we will get independent components (P)
end


% preWhitening code taken from EEGLAB implementation of SOBI by Beloucharani et al.                        (  The paper seems quite exhaustive (P))
% This step makes sure that the data is uncorrelated which is a pre-requisite for ICA. Let us say I have a matrix X and a mixing matrix
% W, and I want to find the matrix W so that cov(X*W) = I.  It turns out that you can find W by finding the matrix equivalent of W= inv(sqrtm(cov(X)));
% This decorrelating process is sometimes called 'whitening' because it produces two uncorrelated signals with unit variance, which is the definition
% of white noise. More info: http://courses.washington.edu/matlab1/matlab/Lesson_17.m
function Q = preWhitenData(Rx)
   IBL=sqrtm(Rx);  % X = sqrtm(A) returns the principal square root of the matrix A, that is, X*X = A. (P)
   Q=inv(IBL);
end


function W_scaled = scaleW(W_unscaled)
    W_scaled = [];  % can be initialized in an efficient way (P)
    for i = 1:1:size(W_unscaled,1)
       W_scaled(i,:) = W_unscaled(i,:)./sqrt(sum(W_unscaled(i,:).^2));  % row wise scaling (P)
    end
end


function printHTML(emg, eog, peaks, pf, ac, el, W, Sps, St, eST, HTMLtitle, cond, freq, toi)           % I didn't read it. MIght ask Chivon if it is required or not
    
    subjid=strrep(HTMLtitle,'.dat','');
    output = strcat(subjid,'smarter_',cond);
    if exist(output,'dir') ~= 7
       fprintf(1, '\nArtifact Source dir doesn''t exist, Creating one...');
       mkdir(output);
    else 
       fprintf(1, '\nArtifact Source dir already exists, Overwriting...'); 
    end
    cd(output);
    A = inv(W);

    % calculating global ERP with all sources intact.
    blockSize = 4096 ;
    Xorig = zeros(size(W,1),blockSize);
    for t = 1:1:size(toi,1)
      %  Xorig = Xorig + A*St(:,toi(t,1):toi(t,2));  
        Xorig = Xorig + A*St(:,(t-1)*blockSize+1:t*blockSize);
    end
    Xorig = Xorig./size(toi,1);
    
    XorigERP = mean(abs(Xorig)); 
    
    
    emg_checked = [];
    
    %determine proper window size for viewing 65hz with current sample rate
    %and NFFT
    increment = freq(2)-freq(1);
    windowfreqxaxis = (120 + increment)/increment; %change 120 to the desired max viewing frequency on x-axis
    windowfreqxaxis = round(windowfreqxaxis);
    
    %open a file that will store the dif for each source%%
    fid = fopen(strcat(strrep(HTMLtitle,'.dat',''),'_SourceDiffs.txt'),'w');
    if fid == -1
        fprintf(1,'Error creating SourceDifs file\n');
    end

    fprintf(fid,'Source\tpercent dif \t\n');

   % printing images
    for i = 1:1:length(emg)
        scrsz = get(0,'ScreenSize'); figure('Visible','off','Position',[1 1 scrsz(3)*.95 scrsz(4)*.15]);
        %figure('Visible','off','Position',[1 1 2432 288]); %creates output based on a 2560 x 1440 screen resolution
        subplot(1,7,1), topoplot(A(:,emg(i)), el, 'electrodes', 'off');
       % subplot(1,7,2), semilogy(Sps(:),freq(1:windowfreqxaxis), Sps(emg(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,7,2),  semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(emg(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);
        %subplot(1,7,3), loglog(freq(1:104), Sps(emg(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,7,3), loglog(freq(1:windowfreqxaxis), Sps(emg(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(emg(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end
        subplot(1,7,4), plot(ac(emg(i),:),'b+'); axis([1 64 -.2 1]); hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        %subplot(1,7,4), plot(ac(emg(i),:),'b+'); axis([4 64 -.2 1]); hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms of the epoch from the autocorrelation display
        if mean(ac(emg(i),:)) < .4
           title(strcat('min1to20lags:',num2str(min(ac(emg(i),1:20)))),'Color','r'); 
%            set(gcf,'Color','c');
           emg_checked = [emg_checked, emg(i)];
        else
           title(strcat('min1to20lags:',num2str(min(ac(emg(i),1:20)))),'Color','b');  
        end
        
        % look at the effect of removing artifact on the evoked response
        % potential (ERP) for all channels.
        blockSize = 4096;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),emg(i)))*St(setdiff(1:size(W,1),emg(i)),(t-1)*blockSize+1:t*blockSize);        
        end
        Xs = Xs./size(toi,1);
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels.        
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight'; 
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        dif = 100*mean(abs(XorigERP - XsERP)./XorigERP); % in percentage increase or decrease;
        diftest = mean((XorigERP - XsERP)./XorigERP); % test to see if it adds up to zero without absolute values calcuations
        
        fprintf(fid, '%d\t%d\t\r\n', emg(i),diftest); %changed to print out new dif calc
  
        if dif >= 1
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','r');        
        else
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','b');        
        end
        
        subplot(1,7,6), plot(eST(emg(i),:));  axis 'tight'; set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        rms_ = sqrt(sum(diff(eST(emg(i),:)).^2)/length(eST(emg(i),:)));
        title(strcat('RMS = ',num2str(rms_))); 
        
        subplot(1,7,7), plot(eST(emg(i),1:2560));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 512; 1536; 2559]); set(gca,'XTickLabel',[-250;0;500;1000])  % 2662:4300
        
                set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('emg_',num2str(emg(i))), 'tiff');
        close all;
        %cropping extra white space.
%         imgs = imread(strcat('emg_',num2str(emg(i)),'.tif')); imwrite(imgs(:,300:end-150,:),strcat('emg_',num2str(emg(i)),'.tif'));
    end
    
    eog_checked = [];
    for i = 1:1:length(eog)
        scrsz = get(0,'ScreenSize'); figure('Visible','off','Position',[1 1 scrsz(3)*.95 scrsz(4)*.15]);
        %figure('Visible','off','Position',[1 1 2432 288]); %creates output based on a 1920 x 1200 screen resolution
        subplot(1,7,1), topoplot(A(:,eog(i)), el, 'electrodes', 'off');
        %subplot(1,7,2), semilogy(freq(1:windowfreqxaxis), Sps(eog(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,7,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(eog(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);       
        
        %subplot(1,7,3), loglog(freq(1:104), Sps(eog(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,7,3), loglog(freq(1:windowfreqxaxis), Sps(eog(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(eog(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end
        subplot(1,7,4), plot(ac(eog(i),:),'b+'); axis([1 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        %subplot(1,7,4), plot(ac(eog(i),:),'b+'); axis([4 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms of the epoch from the autocorrelation display
        if min(ac(eog(i),:)) < .4
           title(strcat('mean:',num2str(mean(ac(eog(i),:)))),'Color','r'); 
        else
           title(strcat('mean:',num2str(mean(ac(eog(i),:)))),'Color','b');  
        end   
        
        % look at the effect of removing artifact on the evoked response
        % potential (ERP) for all channels.
        blockSize = 4096;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),eog(i))) * St(setdiff(1:size(W,1),eog(i)),(t-1)*blockSize+1:t*blockSize);        
        end
        Xs = Xs./size(toi,1);
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels.         
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight';
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        dif = 100*mean(abs(XorigERP - XsERP)./XorigERP); % in percentage increase or decrease;
        diftest = mean((XorigERP - XsERP)./XorigERP); % test to see if it adds up to zero without absolute values calcuations
        fprintf(fid, '%d\t%d\t\r\n', eog(i),diftest);
        
        if mean(ac(eog(i),:)) > 0.9 && (dif >= 1)         %%%%%%greater than 1% change%%%%%%%
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','r'); 
            eog_checked = [eog_checked, eog(i)];
%             set(gcf,'Color','y');  
        elseif mean(ac(eog(i),:)) < 0.9 && (dif >= 1)  
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','w'); 
%             set(gcf,'Color','m');              
        else
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','w');        
%             set(gcf,'Color','b');         %%%%%ignore these%%
        end
        
        subplot(1,7,6), plot(eST(eog(i),:)); axis 'tight'; set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        rms_ = sqrt(sum(diff(eST(eog(i),:)).^2)/length(eST(eog(i),:)));
        title(strcat('RMS = ',num2str(rms_)));   
        
        subplot(1,7,7), plot(eST(eog(i),1:2560));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 512; 1536; 2559]); set(gca,'XTickLabel',[-250;0;500;1000])       %2662:4300 
                set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('eog_',num2str(eog(i))), 'tiff');        
        close all;
        %cropping extra white space.
%         imgs = imread(strcat('eog_',num2str(eog(i)),'.tif')); imwrite(imgs(:,300:end-150,:),strcat('eog_',num2str(eog(i)),'.tif'));
        
    end
    
    peaks_checked = [];
    for i = 1:1:length(peaks)
        scrsz = get(0,'ScreenSize'); figure('Visible','off','Position',[1 1 scrsz(3)*.95 scrsz(4)*.15]);
        %figure('Visible','off','Position',[1 1 2432 288]); %creates output based on a 1920 x 1200 screen resolution
        subplot(1,7,1), topoplot(A(:,peaks(i)), el, 'electrodes', 'off');
        %subplot(1,7,2), semilogy(freq(1:windowfreqxaxis), Sps(peaks(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
       subplot(1,7,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(peaks(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);        
        
        %subplot(1,7,3), loglog(freq(1:104), Sps(peaks(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,7,3), loglog(freq(1:windowfreqxaxis), Sps(peaks(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(peaks(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end        
        
        subplot(1,7,4), plot(ac(peaks(i),:),'b+'); axis([1 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        %subplot(1,7,4), plot(ac(peaks(i),:),'b+'); axis([4 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms of the epoch from the autocorrelation display
        if min(ac(peaks(i),:)) < .4
           title(strcat('min1to20lags:',num2str(min(ac(peaks(i),1:20)))),'Color','r'); 
           peaks_checked = [peaks_checked, peaks(i)];
%            set(gcf,'Color','g');
        else
           title(strcat('min1to20lags:',num2str(min(ac(peaks(i),1:20)))),'Color','b');  
        end        
        
        % look at the effect of removing artifact on the evoked response
        % potential (ERP) for all channels.
        blockSize = 4096;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),peaks(i)))*St(setdiff(1:size(W,1),peaks(i)),(t-1)*blockSize+1:t*blockSize);        
        end
        Xs = Xs./size(toi,1);
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels.        
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight'; set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        dif = 100*mean(abs(XorigERP - XsERP)./XorigERP); % in percentage increase or decrease;
        diftest = mean((XorigERP - XsERP)./XorigERP); % test to see if it adds up to zero without absolute values calcuations
        fprintf(fid, '%d\t%d\t\r\n', peaks(i),diftest);
        
        if dif >= 1
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','r');        
        else
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','b');        
        end

        subplot(1,7,6), plot(eST(peaks(i),:));  axis 'tight'; set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        rms_ = sqrt(sum(diff(eST(peaks(i),:)).^2)/length(eST(peaks(i),:)));
        title(strcat('RMS = ',num2str(rms_)));   
        
        subplot(1,7,7), plot(eST(peaks(i),1:2560)); axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 512; 1536; 2559]); set(gca,'XTickLabel',[-250;0;500;1000])       % 2662:4300
        
        set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('peaks_',num2str(peaks(i))), 'tiff');        
        close all;
        %cropping extra white space.
%         imgs = imread(strcat('peaks_',num2str(peaks(i)),'.tif')); imwrite(imgs(:,300:end-150,:),strcat('peaks_',num2str(peaks(i)),'.tif'));
        
    end    
    normal_checked = [];
    normal = setdiff(1:size(St,1),union(union(emg,peaks),eog));
    for i = 1:1:length(normal)
        scrsz = get(0,'ScreenSize'); figure('Visible','off','Position',[1 1 scrsz(3)*.95 scrsz(4)*.15]);
        %figure('Visible','off','Position',[1 1 2432 288]); %creates output based on a 1920 x 1200 screen resolution
        subplot(1,7,1), topoplot(A(:,normal(i)), el, 'electrodes', 'off');
        %subplot(1,7,2), semilogy(freq(1:windowfreqxaxis), Sps(normal(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,7,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(normal(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);        
        
        %subplot(1,7,3), loglog(freq(1:104), Sps(normal(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,7,3), loglog(freq(1:windowfreqxaxis), Sps(normal(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(normal(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end    
        
        subplot(1,7,4), plot(ac(normal(i),:),'b+'); axis([1 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); 
        %subplot(1,7,4), plot(ac(normal(i),:),'b+'); axis([4 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms from the autocorrelation display
        if min(ac(normal(i),:)) < .4
           title(strcat('min1to20lags:',num2str(min(ac(normal(i),1:20)))),'Color','r'); 
%            set(gcf,'Color','c');
        else
           title(strcat('min1to20lags:',num2str(min(ac(normal(i),1:20)))),'Color','b');  
        end      
        
        
        % look at the effect of removing artifact on the evoked response potential (ERP) for all channels. This is done by re-projectiing the current source signal onto sensor space
        % St contains the output of m sources: m x n_trial*epoch_time
        blockSize = 4096;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),normal(i)))*St(setdiff(1:size(W,1),normal(i)),(t-1)*blockSize+1:t*blockSize);        % projecting the current source into sensor space for the current trial
        end
        Xs = Xs./size(toi,1);  % Size: m x 4096
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels due to this source.   Dim:  1x 4096       
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight';  %XorigERP contains the ERP due to all channels
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        dif = 100*mean(abs(XorigERP - XsERP)./XorigERP); % in percentage increase or decrease;
        diftest = mean((XorigERP - XsERP)./XorigERP); % test to see if it adds up to zero without absolute values calcuations
        fprintf(fid, '%d\t%d\t\r\n', normal(i),diftest);
        
        if dif >= 1
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','r');        
        else
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','b');        
        end
 
        % eST: m x epoch_time
        subplot(1,7,6), plot(eST(normal(i),:)); axis 'tight'; set(gca,'XTick',[1; 512; 4095]); set(gca,'XTickLabel',[-250;0;1750]);
        rms_ = sqrt(sum(diff(eST(normal(i),:)).^2)/length(eST(normal(i),:)));  % rms on the difference of the evoked source signal
        title(strcat('RMS = ',num2str(rms_)));
        
        subplot(1,7,7), plot(eST(normal(i),1:2560));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 512; 1536; 2559]); set(gca,'XTickLabel',[-250;0;500;1000])        %2662:4300
        %         if rms_ > 0.[1
%             set(gcf,'Color','y');   
%             normal_checked = [normal_checked, normal(i)];            
%         end        
        
        set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('normal_',num2str(normal(i))), 'tiff');        
        close all;
        %cropping extra white space.
%         imgs = imread(strcat('normal_',num2str(normal(i)),'.tif')); imwrite(imgs(:,300:end-150,:),strcat('normal_',num2str(normal(i)),'.tif'));
        
    end    
    fclose(fid);
    
   % print the HTML file
    fid = fopen('ArtifactsHTML.html','w');
    fprintf(fid, '\n<HTML>\n<TITLE>%s</TITLE>\n<BODY>\n <form action="http://localhost/cgi-bin/poll.cgi" method="POST">\n<TABLE border=1>\n', strcat(HTMLtitle,'_',cond)); % "http://www.rebol.com/cgi-bin/test-cgi.cgi" works
    fprintf(fid, '<input type="text" name="filename", value=%s>\n',strcat(strrep(HTMLtitle,' ',''),'_',cond));
    fprintf(fid, '<input type="text" name="nsources", value=%s>\n',num2str(size(St,1)));
    
    %store an Excel file that lists the sources in their order of appearance. To be uesd for Voting%%
    sourceorder = [];
    %fileid = fopen(strcat(strrep(HTMLtitle,'.dat',''),'_voteforrecon.xls'),'w');
    
    if fid == -1
        fprintf(1,'Error creating Voting file in Excel\n');
    end
    
    
    emg2 = [emg_checked setdiff(emg,emg_checked)];    
    for i = 1:1:length(emg2)
        if ismember(emg2(i),emg_checked)
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm"><input type="radio" name="%d" value="emg" checked><input type="radio" name="%d" value="Rnorm">', num2str(emg2(i)),emg2(i),emg2(i),emg2(i));
        else
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm" checked><input type="radio" name="%d" value="emg"><input type="radio" name="%d" value="Rnorm">', num2str(emg2(i)),emg2(i),emg2(i),emg2(i));
        end
        fprintf(fid, '\n<IMG SRC=''%s''></TD></TR>', strcat('emg_',num2str(emg2(i)),'.tif')); 
        sourceorder = [sourceorder emg2(i)];
    end

    peaks2 = [peaks_checked setdiff(peaks,peaks_checked)];
    for i = 1:1:length(peaks2)
        if ismember(peaks2(i),peaks_checked)
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm"><input type="radio" name="%d" value="peaks" checked><input type="radio" name="%d" value="Rnorm">', num2str(peaks2(i)),peaks2(i), peaks2(i), peaks2(i));        
        else
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm" checked><input type="radio" name="%d" value="peaks"><input type="radio" name="%d" value="Rnorm">', num2str(peaks2(i)),peaks2(i), peaks2(i), peaks2(i));                    
        end
        fprintf(fid, '\n<IMG SRC=''%s''></TD></TR>', strcat('peaks_',num2str(peaks2(i)),'.tif'));  
%        fprintf(fileid, '%d\r\n',peaks2(i));
        sourceorder = [sourceorder peaks2(i)];
    end
    
    eog2 = [eog_checked setdiff(eog,eog_checked)];
    for i = 1:1:length(eog2)
        if ismember(eog2(i), eog_checked)
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2><input type="radio" name="%d" value="norm"><input type="radio" name="%d" value="eog" checked><input type="radio" name="%d" value="Rnorm"><input type="radio" name="%d" value="Reog">', num2str(eog2(i)), eog2(i), eog2(i), eog2(i), eog2(i));        
        else
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2><input type="radio" name="%d" value="norm" checked><input type="radio" name="%d" value="eog"><input type="radio" name="%d" value="Rnorm"><input type="radio" name="%d" value="Reog">', num2str(eog2(i)), eog2(i), eog2(i), eog2(i), eog2(i));        
        end
        fprintf(fid, '\n<IMG SRC=''%s''></TD></TR>', strcat('eog_',num2str(eog2(i)),'.tif'));  
        sourceorder = [sourceorder eog2(i)];
    end
   
    normal2 = [normal_checked setdiff(normal,normal_checked)];
    for i = 1:1:length(normal2)
        if ismember(normal2(i), normal_checked)
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm"><input type="radio" name="%d" value="eog" checked><input type="radio" name="%d" value="Rnorm">', num2str(normal2(i)), normal2(i), normal2(i), normal2(i));        
        else
            fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm" checked><input type="radio" name="%d" value="emg"><input type="radio" name="%d" value="Rnorm">', num2str(normal2(i)), normal2(i), normal2(i), normal2(i));        
        end
        fprintf(fid, '\n<IMG SRC=''%s''></TD></TR>', strcat('normal_',num2str(normal2(i)),'.tif'));  
        sourceorder = [sourceorder normal2(i)];
    end
    fprintf(fid, '\n</TABLE><input type="submit" value="Vote"></BODY></HTML>');
    
    fclose(fid);
    xlswrite(strcat(strrep(HTMLtitle, '.dat',''),'_voteforrecon.xlsx'),sourceorder');
    cd('..');    
end

