% This code is developed to run SOBI on EEG data; 
% author: Manish Saggar (mishu@cs.utexas.edu). 
% Iman has changed the code and added the embedded source trial  ver:04/12/2012
% (irezazadeh@ucdavis.edu)
%Chivon has changed the code to run with CPT task data. 

function y = RegenerateSourceProjections(eegfilename, eegfile, outputDir)
 
    
    eegfile = strtrim(eegfile);        
    
    % Variables
    hdr = readBESAsb_header([eegfile '_epoched_forSOBI']);
    load([outputDir '_epoched_forSOBI\' eegfilename '_epoched_forSOBI_Sources_CPT\' eegfilename '_epoched_forSOBI_CPT_W.mat']);
          W_scaled = real(scaleW(W));    

            % output nSources x nTrials per condition; S = WX;            
            createSourceTrials(hdr, W_scaled, eegfile, toi, eegfilename, trig);

            % output projection of each source in sensor space for the 
            % averaged epoch;            
           createEvokedSourceTrial(hdr, W_scaled, eegfile, toi, eegfilename);
           createEmbeddedEvokedSourceTrial(hdr, W_scaled, eegfile, toi, eegfilename, trig);
            
end

function sT = createEvokedSourceTrial(hdr, W, eegfile, toi, eegfilename)
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
    writeBESAsb_data(sT, strcat(eegfilename,'_CPT_EvokedSourceTrial.dat'));
    
    % creating event file
    fid = fopen(strcat(eegfilename,'_CPT_EvokedSourceTrial.evt'),'w');
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
    writeBESAsb_data(chData, strcat(eegfilename,'_CPT_EvokedSourceTrialSensorProj.dat'));
 % creating event file
    fid = fopen(strcat(strrep(eegfile,'.dat',''),'_CPT_EvokedSourceTrialSensorProj.evt'),'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Sensor Projection\n');
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    for t = 1:1:size(W,1)
       fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*(blockSize/2048))*1e6+1, 1, t, strcat('Trig.',num2str(t)));
    end
    fclose(fid);   
       
end
function sT = createEmbeddedEvokedSourceTrial(hdr, W, eegfile, toi, evtfile, eegfilename)
    blockSize = 5530;
    no_event=2;
    sT_new = zeros(size(W,1), blockSize,no_event);
    %sT = zeros(size(W,1), blockSize*no_event);
    sT = [];
    trigger_counter=zeros(no_event);
    
    evts1 = textread(evtfile,'%s');
    tmu = [];
    trig = [];
    for i = 5:1:size(evts1,1)
        if mod(i,5) == 0
            tmu(end+1) = str2num(cell2mat(evts1(i)));
        elseif mod(i,5) == 2
            trig(end+1) = str2num(cell2mat(evts1(i)));
        end
    end 
    
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
    writeBESAsb_data(sT, strcat(eegfilename,'_CPT_EmbeddedEvokedSourceTrial.dat'));
    
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
    writeBESAsb_data(chData, strcat(eegfilename,'_CPT_EmbeddedEvokedSourceTrialSensorProj.dat'));
        
    % creating event files: one with conditions and another with sources
    fid = fopen(strcat(eegfilename,'_CPT_EmbeddedEvokedSourceTrialSensorProj.evt'),'w');
    fid2 = fopen(strcat(eegfilename,'_CPT_EmbeddedEvokedSourceTrialSensorProjSourceNums.evt'),'w');
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
function sT = createSourceTrials(hdr, W, eegfile, toi, eegfilename, trig)
    subjid=eegfilename;    
    output = strcat(subjid,'_Sources_CPT');
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
          
           
        end 
           writeBESAsb_data(chData_tmp, strcat(eegfilename,'_CPT_Source_',num2str(s),'_Projection_Into_SensorSpace_AllTrials.dat'));    
           clear chData_tmp
    end 
    % creating event file
    fid = fopen(strcat(eegfilename,'_CPT_SensorSpace_AllTrials.evt'),'w');
    fid2 = fopen(strcat(eegfilename,'_CPT_SensorSpace_AllTrialsTrigs.evt'),'w');
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
    writeBESAsb_data(sT, strcat(eegfilename,'_CPT_SourceTrials.dat'));

    % creating event file
    fid = fopen(strcat(eegfilename,'_CPT_SourceTrials.evt'),'w');
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
function W_scaled = scaleW(W_unscaled)
    W_scaled = [];
    for i = 1:1:size(W_unscaled,1)
       W_scaled(i,:) = W_unscaled(i,:)./sqrt(sum(W_unscaled(i,:).^2));
    end
end
