% This code is developed to run SOBI on EEG data; 
% author: Manish Saggar (mishu@cs.utexas.edu). 
% Iman has changed the code and added the embedded source trial  ver:04/12/2012
% (irezazadeh@ucdavis.edu)
%Chivon has changed the code to run with CPT task data. 

function y = tSOBI_CPT_betterSMART_withdetrend(eegfilename, eegfile, evtfile, sfpfile, outputDir, chunks)
 
    diary tSOBI_logFile;
    eegfile = strtrim(eegfile);        
    
    % Variables
    hdr = readBESAsb_header(eegfile);
    tang_tau = [1:10,12:2:20,25:5:100,120:20:300]; % in msec 

    
    %  derive tau by changing tang_tau according to our Fs
    tau = [0 unique(round(hdr.Fs*(tang_tau./1000)))];

    sinThreshold = 1e-08;
 
    
   % Creating trial of interest (toi)variable by reading data from the Event file.
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
origISI = [];
origBlock = [];
origRT = [];
origCount = [];

for i = 9:9:size(evts,1)
    origtmu(end+1) = str2num(cell2mat(evts(i)));
    origtrig(end+1)= str2num(cell2mat(evts(i+2)));
    origISI(end+1) = str2num(cell2mat(evts(i+6)));
    origBlock(end+1) = str2num(cell2mat(evts(i+5)));
    origRT(end+1) = str2num(cell2mat(evts(i+7)));
    origCount(end + 1) = str2num(cell2mat(evts(i+8)));
end

Trialstep = [];
%read the original triggers into tmu and trig variables and mark where
%missing trials were.
for i = 1:1:size(origtrig,2)
    if origtrig(i) >100 && origtrig(i)<317
        trig(end+1) = 2;
        tmu(end+1) = origtmu(i);
%         if origCount(i) == 1
%         Trialstep(end+1) = 1;
%         elseif origCount(i) == origCount(i-1)+1
%             Trialstep(end+1) = 1;
%         elseif origtrig(i-1) > 1             
%             Trialstep(end+1) = origCount(i) - origCount(i-1);
%         else Trialstep(end+1) = origCount(i) - origCount(i-2);
%         end
    elseif origtrig(i) >500 && origtrig(i)<525
        trig(end+1) = 4;
        tmu(end+1) = origtmu(i);
%          if origCount(i) == 1
%         Trialstep(end+1) = 1;
%         elseif origCount(i) == origCount(i-1)+1
%             Trialstep(end+1) = 1;
%         elseif origtrig(i-1) > 1             
%             Trialstep(end+1) = origCount(i) - origCount(i-1);
%         else Trialstep(end+1) = origCount(i) - origCount(i-2);
%         end
    end
end
     

    toiAll = [];
    trigAll = trig;
    cAll = 1;
    tmu = tmu./1e6; % converting to seconds to multiply by # of samples.
    for i = 1:1:size(trig,2)
        toiAll(cAll,1) = round(tmu(i)*2048) - 3072; %Set the stim-locked epoch range in samples         
        toiAll(cAll,2) = round(tmu(i)*2048) + 2458;
        cAll = cAll + 1;
    end
    
      
    % dividing the toiAll into chunks based upon the value of parameter
    % chunks. This parameter is always set to 1 for the CPT data.
    
  
    if chunks == 1

        % calc correlation matrices
        Rx = createCorrMat(hdr, eegfile, toiAll, tau);

        % Feed Correlation Matrices to Joint Diagonalization.
        [W, B, D] = jointDiagonalization(Rx, sinThreshold);
        W_scaled = real(scaleW(W));    
       
        % output nSources x nTrials per condition; S = WX;            
        
        aST = createSourceTrials(hdr, W_scaled, eegfile, toiAll, 'CPT', eegfilename, trig);

        % output projection of each source in sensor space for the 
        % averaged epoch;            
       
       eST = createEvokedSourceTrial(hdr, W_scaled, eegfile, toiAll, 'CPT', eegfilename);

        embeddedST = createEmbeddedEvokedSourceTrial(hdr, W_scaled, eegfile, toiAll, 'CPT',evtfile,eegfilename,trig,tmu);

        % create HTML file using SMART
           smart_autism(aST, eST, hdr, W_scaled, eegfilename, toiAll, trigAll, 'CPT', sfpfile);
    else
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
function smart_autism(aST, eST, hdr, W, eegfile, toi, trig, cond, sfpfile)
    emgThreshold = 0.6;
    eogThreshold = 0.7;
    el = readlocs(sfpfile);
    % calculate EMG and EOG components based on Auto-correlation function on the
    % evoked (average) source trial.
    [emgnew, ac] = findACFbasedEMG(aST, W, toi, emgThreshold);     
    eognew = findACFbasedEOG(ac, W, eogThreshold);
    % calculate Peak components based on PowerSpectrum analysis on each
    % Source trial. 
    [peaksnew peaks Sps freq]= findPeaks(aST, hdr, W, toi);  
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path
    printHTML(emgnew, eognew, peaksnew, peaks, ac, el, W, Sps, aST, eST,...
        eegfile, cond, freq, toi);    
    save(strcat(strrep(eegfile,'.dat',''),'_',cond,'_W'),'W','aST','toi','trig');
end
function [y acSt] = findACFbasedEMG(St, W, toi, thresh)
    y = [];
    acSt = [];
    blockSize = 5530;
    for s = 1:1:size(W,1)
        acst = [];
        for t = 1:1:size(toi,1)
            %  ac = acf_old(St(s, toi(t,1):toi(t,2))',64);
             ac = acf_old(St(s, (t-1)*blockSize+1:t*blockSize)',64); 
             acst(t,:) = ac.ac;
        end
        % average the acf for all trials for each source.
        acSt(s,:) = mean(acst);
        if min(acSt(s,1:20)) < thresh  % trying min instead of mean
            y = [y s];
        end
    end
end
function [y] = findACFbasedEOG(ac, W, thresh)
    y = [];
    for i = 1:1:size(W,1)
       if mean(ac(i,:)) > thresh
           y = [y i];
       end
    end
end
function [y pp psStNL freq] = findPeaks(St, hdr, W, toi)
    y =[];
    pp = [];
    NFFT = 2700;    % controls the frequency resolution.
    Hs = spectrum.welch({'Hann','periodic'}, NFFT, 0);
    psSt = [];
    psStNL = [];
    blockSize = 5530;
    % calculate power spectrum for each trial and each source.
    for s = 1:1:size(W,1)
        psst = [];
        for t = 1:1:size(toi,1)
           % tmp = St(s, toi(t,1):toi(t,2));
            tmp = St(s, (t-1)*blockSize+1:t*blockSize);
            ps = psd(Hs, tmp,'Fs',hdr.Fs);
            psst(t,:) = ps.Data';
        end
        % average the power spectrum for all trials for each source.
        psSt(s,:) = 10*log10(mean(psst));
        psStNL(s,:) = mean(psst);
        freq = ps.Frequencies;

        % to avoid error from findpeaks
        if max(psStNL(s,:)) <.5, 
            continue;
        end

        [p,l] = findpeaks(psStNL(s,:),'SORTSTR','descend','MINPEAKHEIGHT',.5);
        npeaks = 5;
        if length(l) < 1
            continue;
        elseif length(l) < 5
            npeaks = length(l);
        end

        freqPeaks = round(freq(l));
        for i = 1:1:npeaks   % just look at first three peaks
           if freqPeaks(i) == 59 || freqPeaks(i) == 60 || freqPeaks(i) == 61
                y = [y s];
                pp = [pp freqPeaks(i)];  
           end
        end
    end
end
function sT = createEvokedSourceTrial(hdr, W, eegfile, toi, cond, eegfilename)
    blockSize = 5530;
    sT = zeros(size(W,1), blockSize);
    for t = 1:1:size(toi,1)
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2));
        tmp = W*x;
        % demeaning the sources from -80: 20  msec (this is shifted by 1 sample also).
        tmp = tmp - mean(tmp(:,2908:3113),2)*ones(1,size(tmp,2));
        sT = sT + tmp; 
    end
    sT = sT./size(toi,1);
%     sT = sT * 1e3;
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path    
    writeBESAsb_data(sT, strcat(eegfilename,'_',cond,'_EvokedSourceTrial.dat'));
    
    % creating event file
    fid = fopen(strcat(eegfilename,'_',cond,'_EvokedSourceTrial.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Sensor Projection\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
  %  for t = 1:1:size(W,1)
       fprintf(fid, '%d\t%d\t%d\t%s\n', 1500*1e3, 1, 1, 'Stimulus');
  %  end
    fclose(fid); 
    
    % creating evoked source trial projections
    chData = zeros(size(W,1), size(W,1)*blockSize);
    for t = 1:1:size(W,1)  % for each source
        % append projection of each source into sensor space
        A = inv(W);
        chData(:,(t-1)*blockSize+1:t*blockSize) = A(:,t)*sT(t,:);
    end
    writeBESAsb_data(chData, strcat(eegfilename,'_',cond,'_EvokedSourceTrialSensorProj.dat'));
 % creating event file
    fid = fopen(strcat(strrep(eegfile,'.dat',''),'_',cond,'_EvokedSourceTrialSensorProj.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Sensor Projection\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    for t = 1:1:size(W,1)
       fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048))*1e6+1, 1, t, strcat('Trig.',num2str(t)));
    end
    fclose(fid);   
       
end
function sT = createEmbeddedEvokedSourceTrial(hdr, W, eegfile, toi, cond,evtfile, eegfilename,trig,tmu)
    blockSize = 5530;
    no_event=2;
    sT_new = zeros(size(W,1), blockSize,no_event);
    %sT = zeros(size(W,1), blockSize*no_event);
    sT = [];
    trigger_counter=zeros(no_event);
    
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
        trigger_no=trig(t)/2;
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2));
        tmp = W*x;
        % demeaning the sources from -80: 20  msec (this is shifted by 1 sample also).
        tmp = tmp - mean(tmp(:,2908:3113),2)*ones(1,size(tmp,2));
        sT_new(:,:,trigger_no) = sT_new(:,:,trigger_no) + tmp; 
        trigger_counter(trigger_no)=trigger_counter(trigger_no)+1;
    end
    %sT = sT./size(toi,1);
    for i=1:no_event
        sT_new(:,:,i) = sT_new(:,:,i)./trigger_counter(i);
        sT=[sT sT_new(:,:,i)];
    end
    
%     sT = sT * 1e3;
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path    
    writeBESAsb_data(sT, strcat(eegfilename,'_',cond,'_EmbeddedEvokedSourceTrial.dat'));
    
    % creating evoked source trial projections
    
%    pData = zeros(size(W,1), size(W,1)*blockSize,event_no);
    chData = zeros(size(W,1), size(W,1)*blockSize*no_event);    
    bias_time=0;   
    for t = 1:1:size(W,1)  % for each source
        for i=0:no_event-1
        % append projection of each source into sensor space
        A = inv(W);
        %chData(:,(t+i-1)*blockSize+1+bias_time:(t+i)*blockSize+bias_time) = A(:,t)*sT_new(t,:,i+1);
        chData(:,(i)*blockSize+1+bias_time:(i+1)*blockSize+bias_time) = A(:,t)*sT_new(t,:,i+1);
        end
        bias_time=((1+i)*blockSize+bias_time);
   end 
    writeBESAsb_data(chData, strcat(eegfilename,'_',cond,'_EmbeddedEvokedSourceTrialSensorProj.dat'));
        
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
             fprintf(fid2, '%d\t%d\t%d\t%s\n', ((i)*(blockSize/2048)*1e6+bias_time)+1, 1, t, strcat('Trig.',num2str(t)));
            fprintf(fid, '%d\t%d\t%d\t%s\n', ((i)*(blockSize/2048)*1e6+bias_time)+1, 1, (i+1)*2, strcat('Trig.',num2str(t)));
       end 
       bias_time=((i+1)*(blockSize/2048))*1e6+bias_time;
    end
%     for t = 1:1:size(trig,2)
%        fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048))*1e6+1, 1, trig(t), strcat('Trial:',num2str(t)));
%     end
    fclose(fid);
    fclose(fid2);    
end 
function sT = createSourceTrials(hdr, W, eegfile, toi, cond, eegfilename, trig)
    subjid=eegfilename;    
    output = strcat(subjid,'_Sources_',cond);
    if exist(output,'dir') ~= 7
       fprintf(1, '\nArtifacts dir doesn''t exist, Creating one...');
       mkdir(output);
    else 
       fprintf(1, '\nArtifacts dir already exists, Overwriting...'); 
       mkdir(output);
    end
    cd(output);    
    blockSize = toi(1,2)- toi(1,1);
    sT = zeros(size(W,1), size(toi,1)*(toi(1,2)- toi(1,1)));
    for t = 1:1:size(toi,1)
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2));
        tmp = W*x;
        % demeaning the sources from -80: 20 msec (shifted by 1 sample).
        tmp = tmp - mean(tmp(:,2908:3113),2)*ones(1,size(tmp,2));
        %detrending the sources by epoch.
         tmp = detrend(tmp,'linear');
        sT(:, (t-1)*blockSize+1:t*blockSize) = tmp;
    end
    A = inv(W);
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
     
    eegfile = tmp2; % without full path  
    clear tmp tmp1 tmp2 tmp3;
    for s = 1:1:size(W,1)  % for each source
        chData_tmp = zeros(size(W,1), size(toi,1)*blockSize);    
        for t = 1:1:size(toi,1)
              
           chData_tmp(:,(t-1)*blockSize+1:t*blockSize)= A(:,s)*sT(s,(t-1)*blockSize+1:t*blockSize);
           %chData_tmp(:,(t-1)*blockSize+1:(t*blockSize)-(blockSize/2))= A(:,s)*sT(s,(t-1)*blockSize+1:(t*blockSize)-(blockSize/2));
           %chData_tmp(:,(t*blockSize)-(blockSize/2)+1:t*blockSize)= A(:,s)*sT(s,(t*blockSize)-(blockSize/2)+1:t*blockSize);
           
        end 
           writeBESAsb_data(chData_tmp, strcat(eegfilename,'_',cond,'_Source_',num2str(s),'_Projection_Into_SensorSpace_AllTrials.dat'));    
           clear chData_tmp
    end 
    % creating event file
    fid = fopen(strcat(eegfilename,'_',cond,'_SensorSpace_AllTrials.evt'),'w');
    fid2 = fopen(strcat(eegfilename,'_',cond,'_SensorSpace_AllTrialsTrigs.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Source Trials\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    fprintf(fid2,'Tmu\tCode\tTriNo\tComnt\n');
    for t = 1:1:size(toi,1)
         %fprintf(fid, '%d\t%d\t%d\t%s\n', toi(t,1)*1e3, 1, t, strcat('Trial:',num2str(t)));
        fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048)*1e6+1), 1, t, strcat('Trial:',num2str(t)));
        fprintf(fid2, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048)*1e6+1), 1, trig(t), strcat('Trial:',num2str(t)));
    end
    fclose(fid);
    fclose(fid2);
%     sT = sT * 1e3;    
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path    
    writeBESAsb_data(sT, strcat(eegfilename,'_',cond,'_SourceTrials.dat'));

    % creating event file
    fid = fopen(strcat(eegfilename,'_',cond,'_SourceTrials.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Source Trials\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    for t = 1:1:size(toi,1)
         %fprintf(fid, '%d\t%d\t%d\t%s\n', toi(t,1)*1e3, 1, t, strcat('Trial:',num2str(t)));
        fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048)*1e6+1), 1, t, strcat('Trial:',num2str(t)));
    end
    fclose(fid);
end
function writeBESAsb_data(X_cap, new_file)
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
function Rx = createCorrMat(hdr, eegfile, toi, tau)
    Rx = zeros(hdr.nChans, hdr.nChans, length(tau));
    %blockSize = 2999;
    % Read data
    for t = 1:1:size(toi,1)
        % the data should be read in nChans x Time format.
        %x = readBESAsb_data(eegfile, hdr, toi(t,1)*2+1, toi(t,2)*2 + max(tau));
        x = readBESAsb_data(eegfile, hdr, toi(t,1), toi(t,2) + max(tau));
        % de-mean data
        if isnan(x)~=1
        x = x - mean(x,2)*ones(1, size(x,2));
        
               % calculating corr matrix at different tau values
        for i = 1:1:length(tau)
          %  tt=size(x); 
          %  if tt(1,2)==toi 
                
               Rxx = x(:,1:toi(t,2)-toi(t,1)) * x(:,1+tau(i):toi(t,2)- toi(t,1)+tau(i))';
              % Rxx = x(:,1:blockSize) * x(:,1+tau(i):blockSize+tau(i))'; 
               Rx(:,:,i) = Rx(:,:,i) + 0.5*(Rxx + Rxx'); 
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
    Rx = Rx./size(toi,1);    
end
function [W, B, D] = jointDiagonalization(Rx, sinThreshold)
    
    % calculating pre-sphering matrix
    B = preWhitenData(Rx(:,:,1));
    
    % whitening the correlation matrices
    Rx_white = zeros(size(Rx,1),size(Rx,2),size(Rx,3)-1);  % leaving the instantaneuous correlation.
    for i = 1:1:size(Rx,3)-1
      Rx_white(:,:,i) = B * Rx(:,:,i+1) * B';
    end    
    
    % calculating joint diagolization matrix
    fprintf(1,'\n====Joint Diagonalization of Correlation Matrices====\n');
    [W, D] = joint_diag(reshape(Rx_white,size(Rx_white,1),size(Rx_white,2)*size(Rx_white,3)), sinThreshold);
    W = W' * B;
end
% preWhitening code taken from EEGLAB implementation of SOBI by
% Beloucharani et al.
function Q = preWhitenData(Rx)
   IBL=sqrtm(Rx);
   Q=inv(IBL);
end
function W_scaled = scaleW(W_unscaled)
    W_scaled = [];
    for i = 1:1:size(W_unscaled,1)
       W_scaled(i,:) = W_unscaled(i,:)./sqrt(sum(W_unscaled(i,:).^2));
    end
end
function printHTML(emg, eog, peaks, pf, ac, el, W, Sps, St, eST, HTMLtitle, cond, freq, toi)
    
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
    blockSize = 5530;
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
        subplot(1,7,2),  semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(emg(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; 60; 120; 241]); set(gca,'XTickLabel',[1;30;60;120]);
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
        blockSize = 5530;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),emg(i)))*St(setdiff(1:size(W,1),emg(i)),(t-1)*blockSize+1:t*blockSize);        
        end
        Xs = Xs./size(toi,1);
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels.        
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight'; 
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
        dif = 100*mean(abs(XorigERP - XsERP)./XorigERP); % in percentage increase or decrease;
        diftest = mean((XorigERP - XsERP)./XorigERP); % test to see if it adds up to zero without absolute values calcuations
        
        fprintf(fid, '%d\t%d\t\r\n', emg(i),diftest); %changed to print out new dif calc
  
        if dif >= 1
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','r');        
        else
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','b');        
        end
        
        subplot(1,7,6), plot(eST(emg(i),:));  axis 'tight'; set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
        rms_ = sqrt(sum(diff(eST(emg(i),:)).^2)/length(eST(emg(i),:)));
        title(strcat('RMS = ',num2str(rms_))); 
        
        subplot(1,7,7), plot(eST(emg(i),2662:4300));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 410;1024;1638]); set(gca,'XTickLabel',[-200;0;300;600]);
        
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
        subplot(1,7,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(eog(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; 60; 120; 241]); set(gca,'XTickLabel',[1;30;60;120]);        
        
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
        blockSize = 5530;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),eog(i)))*St(setdiff(1:size(W,1),eog(i)),(t-1)*blockSize+1:t*blockSize);        
        end
        Xs = Xs./size(toi,1);
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels.         
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight';
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
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
        
        subplot(1,7,6), plot(eST(eog(i),:)); axis 'tight'; set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
        rms_ = sqrt(sum(diff(eST(eog(i),:)).^2)/length(eST(eog(i),:)));
        title(strcat('RMS = ',num2str(rms_)));   
        
        subplot(1,7,7), plot(eST(eog(i),2662:4300));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 410;1024;1638]); set(gca,'XTickLabel',[-200;0;300;600])        
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
       subplot(1,7,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(peaks(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; 60; 120; 241]); set(gca,'XTickLabel',[1;30;60;120]);        
        
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
        blockSize = 5530;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),peaks(i)))*St(setdiff(1:size(W,1),peaks(i)),(t-1)*blockSize+1:t*blockSize);        
        end
        Xs = Xs./size(toi,1);
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels.        
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight'; set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        dif = 100*mean(abs(XorigERP - XsERP)./XorigERP); % in percentage increase or decrease;
        diftest = mean((XorigERP - XsERP)./XorigERP); % test to see if it adds up to zero without absolute values calcuations
        fprintf(fid, '%d\t%d\t\r\n', peaks(i),diftest);
        
        if dif >= 1
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','r');        
        else
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','b');        
        end

        subplot(1,7,6), plot(eST(peaks(i),:));  axis 'tight'; set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
        rms_ = sqrt(sum(diff(eST(peaks(i),:)).^2)/length(eST(peaks(i),:)));
        title(strcat('RMS = ',num2str(rms_)));   
        
        subplot(1,7,7), plot(eST(peaks(i),2662:4300));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 410;1024;1638]); set(gca,'XTickLabel',[-200;0;300;600])        
        
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
        subplot(1,7,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(normal(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; 60; 120; 241]); set(gca,'XTickLabel',[1;30;60;120]);        
        
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
        
        
        % look at the effect of removing artifact on the evoked response
        % potential (ERP) for all channels.
        blockSize = 5530;
        Xs = zeros(size(W,1),blockSize);
        for t = 1:1:size(toi,1)
            Xs = Xs + A(:,setdiff(1:size(W,1),normal(i)))*St(setdiff(1:size(W,1),normal(i)),(t-1)*blockSize+1:t*blockSize);        
        end
        Xs = Xs./size(toi,1);
        XsERP = mean(abs(Xs));  % calculating global ERP on all channels.            
                
        subplot(1,7,5), plot(XorigERP, 'b'); hold on; plot(XsERP,'r'); axis 'tight';
        set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
        dif = 100*mean(abs(XorigERP - XsERP)./XorigERP); % in percentage increase or decrease;
        diftest = mean((XorigERP - XsERP)./XorigERP); % test to see if it adds up to zero without absolute values calcuations
        fprintf(fid, '%d\t%d\t\r\n', normal(i),diftest);
        
        if dif >= 1
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','r');        
        else
            title(strcat('R(X),B(Xo),Dif(%)=',num2str(dif)),'Color','b');        
        end
 
        subplot(1,7,6), plot(eST(normal(i),:)); axis 'tight'; set(gca,'XTick',[1; 3072; 5529]); set(gca,'XTickLabel',[-1500;0;1200]);
        rms_ = sqrt(sum(diff(eST(normal(i),:)).^2)/length(eST(normal(i),:)));
        title(strcat('RMS = ',num2str(rms_)));   
subplot(1,7,7), plot(eST(normal(i),2662:4300));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 410;1024;1638]); set(gca,'XTickLabel',[-200;0;300;600])        
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
    fprintf(fid, '\n<HTML>\n<TITLE>%s</TITLE>\n<BODY>\n <form action="http://localhost/cgi-bin/poll.cgi" method="POST">\n<TABLE border=1>\n', strcat(HTMLtitle,'_',cond));
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

