% This code is developed to run SOBI on EEG data; 
% author: Manish Saggar (mishu@cs.utexas.edu). 
% Iman has changed the code and added the embedded source trial  ver:04/12/2012
% (irezazadeh@ucdavis.edu)
%Chivon has changed the code to run with STR task data. 

% Last modified before Pranshul: 18th Feb 2015


function y = tSOBI_CMP_betterSMART(eegfilename, eegfile, evtfile, sfpfile, outputDir, chunks)
 
    diary tSOBI_logFile;   % writes everything printed on the command window
    eegfile = strtrim(eegfile);  % removes leading and trailing whitespace characters from txt       
    
    % Variables
    hdr = readBESAsb_header(eegfile); % this returns a structure having all the information which a generic file has (P)
    tang_tau = [1:10,12:2:20,25:5:100,120:20:300]; % in msec                These are the times for which we will check correlation (P)

    
    % derive tau by changing tang_tau according to our Fs
    tau = [0 unique(round(hdr.Fs*(tang_tau./1000)))]; % this is like the number of sample points for the time points mentioned in tang_tau (P)

    sinThreshold = 1e-08;         % Threshold for non-diagonal extent during SOBI (P)
 
    
   % Creating trial of interest (toi) variable by reading data from the Event file.
    fid = fopen(evtfile,'r');
    MyText = textscan(fid,'%d %d %d','headerlines',1);
    fclose(fid);

    origtmu = MyText{1,1}(:); % includes artifacts' timestamp
    origtrig = MyText{1,3}(:);  % includes arifacts
    tmu = [];
    trig = [];

    start = 1;
    while(origtrig(start) ~= 13)   % 13 is the code for staring of the introduction of compassion meditation
        start = start +1;
    end

    toiAll = [origtmu(start), origtmu(end)];
    j = 0;
    artifact_interval = origtmu(start);

    for i = start: size(origtmu,1)

        if (origtrig(i) == 21 && origtrig(i+1) == 22)
            toiAll(end,2) = origtmu(i);
            toiAll(end+1,:) = [origtmu(i+1),origtmu(end)];
            artifact_interval = artifact_interval + (origtmu(i+1) - origtmu(i));      
        elseif (origtrig(i) == 21 && origtrig(i+1) ~= 22)
            for j = i:size(origtmu,1)
                if origtrig(j) == 22          % this j holds the index for artifact end
                    break;
                end
            end
            artifact_interval_old = artifact_interval;
            artifact_interval = artifact_interval + (origtmu(j) - origtmu(i)); 
            for k = i+1:j-1
                if any(origtrig(k) == [10,13,14,15,16,17,18,19])  % for start conditions
                    trig(end+1) = origtrig(k);   
                    tmu(end+1) = origtmu(j) - artifact_interval; % they will now start at the end of the artifact
                elseif any(origtrig(k) == [20,23,24,25,26,27,28,29])  % for end condition
                    trig(end+1) = origtrig(k);   
                    tmu(end+1) = origtmu(i) - artifact_interval_old; % they will now start at the beginning of the artifact
                else
                    fprintf('Invalid trigger')
                end
            end

            toiAll(end,2) = origtmu(i);
            toiAll(end+1,:) = [origtmu(j),origtmu(end)];

        elseif (i>j) && (origtrig(i) ~= 22)  % we don't want the trigger 22 and triggers between arifacts to get considered again
            trig(end+1) = origtrig(i);
            tmu(end+1) = origtmu(i) - artifact_interval + 100; % 100 is added so that it does not start from absolute 0
        end
    end
    
    toiAll = round((2048/1e6)*toiAll);  % conveting it into the sampling points
    %toiAll(:,2) = toiAll(:,2) - 650;   % covering up for tau
    
    toi_diff = toiAll(:,2) - toiAll(:,1);

    toi_stitch = zeros(size(toiAll));    % it has everything in continuity
    toi_stitch(1,1) = 1;
    toi_stitch(1,2) = toi_stitch(1,1) + toi_diff(1);
    for i = 2:size(toiAll,1)
        toi_stitch(i,1) = toi_stitch(i-1,2) + 1;
        toi_stitch(i,2) = toi_stitch(i,1) + toi_diff(i);
    end

    
    tmu = tmu';
    trig = trig';
    trigAll = trig;
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
        aST = createSourceTrials(hdr, W_scaled, eegfile, toiAll, 'CMP', eegfilename, trig); %m x n_clean_time_points

        % averaged source components across all epoch;
        eST = aST; % since we don't have any trials here
        
        % create HTML file using SMART
        smart_autism(aST, eST, hdr, W_scaled, eegfilename, toiAll, trigAll, 'CMP', sfpfile, outputDir, toi_stitch);
        % The above command writes aST, eST, new EOG and EMG artifact information, new peaks around power noise,... all that in the HTML file (P)
        
    end
    
    
    % writing the integrated event file
    filename =strcat('recon_',eegfilename,'_withRespTrigs.evt');
    fid = fopen(filename,'w');
    
    if fid == -1
        fprintf(1,'Error creating event file with integrated response triggers.\n');
    end
    
    fprintf(fid,'Tmu\tCode\tTriNo\n');
    
    for j = 1:size(trig,1)
        fprintf(fid, '%d\t%d\t%d\n', tmu(j), 1, trig(j));
    end
    
    fclose(fid);
    
    
    
    diary off;             
end


function smart_autism(aST, eST, hdr, W, eegfile, toi, trig, cond, sfpfile, outputDir, toi_stitch)  % [aST, eST, hdr, W_scaled, eegfilename, toiAll, trigAll, 'STR', sfpfile]
    emgThreshold = 0.6;   % How are these thresholds decided? Is this what Cliff was talking about as 'classifier by Manish'?
    eogThreshold = 0.7;   % How are these thresholds decided? Is this what Cliff was talking about as 'classifier by Manish'?
   
    update_sfp(sfpfile);   % removes the unnecessary rows in the beginning
    
    [tmp1, tmp2, tmp3] = fileparts(sfpfile); 
    newsfpfile = strcat(tmp1,'\',tmp2,'_updated',tmp3);
   
    el = readlocs(newsfpfile);  % reads electrode location coordinates and other information from a file. (P)
    %output is a structure containing the channel names and locations (if present). It has three fields: 'eloc.labels', 'eloc.theta' and 'eloc.radius' (P)
    % I get this warning: coordinate conversion failed. However, the output seems fine. It has 5 columns: [electrode number, labels, X, Y, Z] (P)
    
    % calculate EMG components based on Auto-correlation function on the evoked (average) source trial.                ACF means autocorrelation function (P)

    [emgnew, ac] = findACFbasedEMG(aST, W, toi, emgThreshold, toi_stitch);   % ac (dim: m X 64) contains mean autocorrelation coeff, and emgnew contains the list of found emg channels (P)
    
    
    % calculate EOG components based on Auto-correlation function on the evoked (average) source trial.          Electrooculography (EOG) is a
    % technique for measuring the corneo-retinal standing potential that exists between the front and the back of the human eye (P)
    eognew = findACFbasedEOG(ac, W, eogThreshold);  % eog new contains the list of found EOG sources (P)
    
    % calculate Peak components based on PowerSpectrum analysis on each Source trial.                  This is kind of finding the bad peaks for power noise around 60 Hz (P)
    [peaksnew peaks Sps freq]= findPeaks(aST, hdr, W, toi, toi_stitch);
    % peaksnew will have the sources indexes of the peaks (P)
    % peaks will have the values of the frequencies at which peak happens (P)
    % Sps will have the mean power spectrum of all sources. Dim: 88 x nFFT (P)
    % freq will have the list of frequencies for which we have the spectra (P)
    
    [tmp1, tmp2, tmp3] = fileparts(eegfile); 
    eegfile = tmp2; % without full path
    printHTML(emgnew, eognew, peaksnew, peaks, ac, el, W, Sps, aST, eST,...
        eegfile, cond, freq, toi, outputDir); 
    
    save(strcat(strrep(eegfile,'.dat',''),'_',cond,'_W'),'W','aST','toi','trig');   % It will create the file and store 'W','aST','toi','trig' in that (P)
end


function [y acSt] = findACFbasedEMG(St, W, toi, thresh, toi_stitch)  % ACF means autocorrelation function (P)
    y = [];
    acSt = [];

    for s = 1:1:size(W,1)  % size(W,1) = m
        acst = []; 
        
        for t = 1:1:size(toi,1)  % n_epochs
            %  ac = acf_old(St(s, toi(t,1):toi(t,2))',64);
             ac = acf_old(St(s, toi_stitch(t,1):toi_stitch(t,2))',64);  % Estimate the coefficients of the autocorrelation (covariance with its own lagged value). 64 is the max lag considered in digital values (P)
             acst(t,:) = ac.ac; % autocorrelation coefficients. Final size of acst would be: n_epochs x 64   (P)
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


function [y pp psStNL freq] = findPeaks(St, hdr, W, toi, toi_stitch)
    y =[];
    pp = [];
    NFFT = 2700;    % controls the frequency resolution.             Thi is being used as segment length, not NFFT (P)        
    % NFFT specifies the number of discrete Fourier transform (DFT) points to use in the PSD estimate. The default nfft is the greater of 256 or the
    % next power of 2 greater than the length of the segments (P). 
    
    %Hs = spectrum.welch({'Hann','periodic'}, NFFT, 0);  % NFFT act as segment length here. 0 overlap (P)
    psSt = [];
    psStNL = [];
    
    % calculate power spectrum for each trial and each source.
    for s = 1:1:size(W,1)   % 88 (P)
        psst = [];
        for t = 1:1:size(toi,1)  % n_trials = 960
           % tmp = St(s, toi(t,1):toi(t,2));
            tmp = St(s, toi_stitch(t,1):toi_stitch(t,2));
            %ps = psd(Hs, tmp,'Fs',hdr.Fs);          % I will have to use pwelch() because psd() is deprecated. ps should be of dim: (2700/2 + 1) x 1 = 1351 x 1
            
            % I will have to use one of the following commands. Both gives the same answer. NFFT is being used as window length as well as the NFFT (P)
            %[pxx, f] = pwelch(a,hann(NFFT,'periodic'),0, NFFT,2048 ); (P)
            %[toi_stitch(t,1),toi_stitch(t,2)]
            %size(tmp)
            n_seg = size(tmp,2);
            [pxx, f] = pwelch(tmp,hanning(n_seg,'periodic'),0, NFFT,2048 ); %(P)
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
 
    sT = [];
    
    for t = 1:1:size(toi,1)
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1), toi(t,2)); % 88 x n_points_for_1_epoch (P)
        tmp = W*x;    % independent source signal. 88 x n_points. So, SOBI is done separately for each segment (P)
        % demeaning the sources from -80: 20 msec (shifted by 1 sample).
        sT = [sT tmp];  % storing independent source signal for all epochs (P) 
    end

    % writing the source space data into files
    %     sT = sT * 1e3;    
    [tmp1, tmp2, tmp3] = fileparts(eegfile);
    eegfile = tmp2; % without full path    
    writeBESAsb_data(sT, strcat(eegfilename,'_',cond,'_SourceTrials.dat'));   %writes .dat and .generic files for the source_space data


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
    toi(:,2) = toi(:,2) - 650; % this is done to accommodate tau
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


function printHTML(emg, eog, peaks, pf, ac, el, W, Sps, St, eST, HTMLtitle, cond, freq, toi,outputDir)           % I didn't read it. MIght ask Chivon if it is required or not
    
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
        subplot(1,4,1), topoplot(A(:,emg(i)), el, 'electrodes', 'off');
       % subplot(1,7,2), semilogy(Sps(:),freq(1:windowfreqxaxis), Sps(emg(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,4,2),  semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(emg(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);
        %subplot(1,7,3), loglog(freq(1:104), Sps(emg(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,4,3), loglog(freq(1:windowfreqxaxis), Sps(emg(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(emg(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end
        subplot(1,4,4), plot(ac(emg(i),:),'b+'); axis([1 64 -.2 1]); hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        %subplot(1,7,4), plot(ac(emg(i),:),'b+'); axis([4 64 -.2 1]); hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms of the epoch from the autocorrelation display
        if mean(ac(emg(i),:)) < .4
           title(strcat('min1to20lags:',num2str(min(ac(emg(i),1:20)))),'Color','r'); 
%            set(gcf,'Color','c');
           emg_checked = [emg_checked, emg(i)];
        else
           title(strcat('min1to20lags:',num2str(min(ac(emg(i),1:20)))),'Color','b');  
        end
            
        set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('emg_',num2str(emg(i))), 'tiff');
        close all;       
    end
    
    eog_checked = [];
    for i = 1:1:length(eog)
        scrsz = get(0,'ScreenSize'); figure('Visible','off','Position',[1 1 scrsz(3)*.95 scrsz(4)*.15]);
        %figure('Visible','off','Position',[1 1 2432 288]); %creates output based on a 1920 x 1200 screen resolution
        subplot(1,4,1), topoplot(A(:,eog(i)), el, 'electrodes', 'off');
        %subplot(1,7,2), semilogy(freq(1:windowfreqxaxis), Sps(eog(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,4,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(eog(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);       
        
        %subplot(1,7,3), loglog(freq(1:104), Sps(eog(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,4,3), loglog(freq(1:windowfreqxaxis), Sps(eog(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(eog(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end
        subplot(1,4,4), plot(ac(eog(i),:),'b+'); axis([1 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        %subplot(1,7,4), plot(ac(eog(i),:),'b+'); axis([4 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms of the epoch from the autocorrelation display
        if min(ac(eog(i),:)) < .4
           title(strcat('mean:',num2str(mean(ac(eog(i),:)))),'Color','r'); 
        else
           title(strcat('mean:',num2str(mean(ac(eog(i),:)))),'Color','b');  
        end   
       
        subplot(1,7,7), plot(eST(eog(i),1:2560));  axis 'tight'; set(gca,'XMinorTick','off'); set(gca,'XTick',[1; 512; 1536; 2559]); set(gca,'XTickLabel',[-250;0;500;1000])       %2662:4300 
        set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('eog_',num2str(eog(i))), 'tiff');        
        close all;
    end
    
    peaks_checked = [];
    for i = 1:1:length(peaks)
        scrsz = get(0,'ScreenSize'); figure('Visible','off','Position',[1 1 scrsz(3)*.95 scrsz(4)*.15]);
        %figure('Visible','off','Position',[1 1 2432 288]); %creates output based on a 1920 x 1200 screen resolution
        subplot(1,4,1), topoplot(A(:,peaks(i)), el, 'electrodes', 'off');
        %subplot(1,7,2), semilogy(freq(1:windowfreqxaxis), Sps(peaks(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
       subplot(1,4,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(peaks(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);        
        
        %subplot(1,7,3), loglog(freq(1:104), Sps(peaks(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,4,3), loglog(freq(1:windowfreqxaxis), Sps(peaks(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(peaks(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end        
        
        subplot(1,4,4), plot(ac(peaks(i),:),'b+'); axis([1 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        %subplot(1,7,4), plot(ac(peaks(i),:),'b+'); axis([4 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms of the epoch from the autocorrelation display
        if min(ac(peaks(i),:)) < .4
           title(strcat('min1to20lags:',num2str(min(ac(peaks(i),1:20)))),'Color','r'); 
           peaks_checked = [peaks_checked, peaks(i)];
%            set(gcf,'Color','g');
        else
           title(strcat('min1to20lags:',num2str(min(ac(peaks(i),1:20)))),'Color','b');  
        end        
        set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('peaks_',num2str(peaks(i))), 'tiff');        
        close all;
        
    end    
    
    normal_checked = [];
    normal = setdiff(1:size(St,1),union(union(emg,peaks),eog));
    for i = 1:1:length(normal)
        scrsz = get(0,'ScreenSize'); figure('Visible','off','Position',[1 1 scrsz(3)*.95 scrsz(4)*.15]);
        %figure('Visible','off','Position',[1 1 2432 288]); %creates output based on a 1920 x 1200 screen resolution
        subplot(1,4,1), topoplot(A(:,normal(i)), el, 'electrodes', 'off');
        %subplot(1,7,2), semilogy(freq(1:windowfreqxaxis), Sps(normal(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,4,2), semilogy(mean(Sps(:,1:windowfreqxaxis))); hold 'on'; plot(Sps(normal(i),1:windowfreqxaxis),'r'); axis 'tight'; set(gca,'XTick',[1; round(windowfreqxaxis/4); round(windowfreqxaxis/2); round(windowfreqxaxis)]); set(gca,'XTickLabel',[1;30;60;120]);        
        
        %subplot(1,7,3), loglog(freq(1:104), Sps(normal(i),1:104)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        subplot(1,4,3), loglog(freq(1:windowfreqxaxis), Sps(normal(i),1:windowfreqxaxis)); axis 'tight'; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on');
        cof = polyfit(log(freq(18:104)'), log(Sps(normal(i),18:104)), 1); 
        if cof(1) > 0
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','r');
        else
            title(strcat('LLLC16-100Hz:',num2str(cof(1))),'Color','b');
        end    
        
        subplot(1,4,4), plot(ac(normal(i),:),'b+'); axis([1 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); 
        %subplot(1,7,4), plot(ac(normal(i),:),'b+'); axis([4 64 -.2 1]);  hold on; set(gca,'YMinorTick','on'); set(gca,'XMinorTick','on'); %changed to omit the first 337.5ms from the autocorrelation display
        if min(ac(normal(i),:)) < .4
           title(strcat('min1to20lags:',num2str(min(ac(normal(i),1:20)))),'Color','r'); 
%            set(gcf,'Color','c');
        else
           title(strcat('min1to20lags:',num2str(min(ac(normal(i),1:20)))),'Color','b');  
        end      
        
        set(gcf,'PaperPositionMode','auto'); set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat('normal_',num2str(normal(i))), 'tiff');        
        close all;

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
    xlswrite(strcat(outputDir, '\',strrep(HTMLtitle, '.dat',''),'_voteforrecon.xlsx'),sourceorder');
    csvwrite(strcat(strrep(HTMLtitle, '.dat',''),'_voteforrecon.xlsx'),sourceorder');
    cd('..');    
end

