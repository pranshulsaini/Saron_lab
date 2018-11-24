% This code is written to do reconstruction of EEG channel data after
% removing sources based on OUT file generated by the WebServer.
%
% Input to this function is the XLS file containing one column:
% contains name of the OUT file that needs to be reconstructed.
%
% Expectations - (1) Subject's directory is in
% /Users/Shared/autism/SobiGeneratedSubjectFolder/ and (2) OUT files are
% copied from /Library/Webserver/Documents/ to
% /Users/Shared/autism/OutFilesFromSMART/
function reconEEGdata(xlsfilename, reconType)
    % read the XLS file
    [n, t, r] = xlsread(xlsfilename);
    
    % location of subject dirs
    subloc = '/Users/Shared/autism/SobiGeneratedSubjectFolder/';
    cd(subloc);
    
    % location of sfp files
    outloc = '/Users/Shared/autism/OutFilesFromSMART/';
    
    % start from i = 2, since the first row is header row.
    for i = 2:1:size(r,1)
        % go inside each subject's folder
        subj = r{i,1}(1:end-9);  %strsplit('_',r{i,1});
        subjDir = strcat('/Users/Shared/autism/SobiGeneratedSubjectFolder/',subj);
        if exist(subjDir,'dir') ~= 7
           fprintf(1, '\nDirectory for the subject (%s) doesn''t exists, Run tSOBI first',subjDir);
           return;
        else 
           fprintf(1, '\nSubject''s dir already exists, Going in...'); 
        end
        cd(subjDir);  
        wfile = strrep(r{i,1},'.out','_W.mat');
        wf = load(wfile);
        W = wf.W;
        St = wf.aST;
        toi = wf.toi;
        trig = wf.trig;
        reconDatUsingOut(W, St, toi, trig, strcat(outloc,r{i,1}), reconType);     
    end    
end

% Output: reconstructed EEG data.
function reconDatUsingOut(W, St, toi, trig, outfile, reconType)
    A = inv(W);
    fprintf(1,'\n File under consideration = %s\n',outfile);
    if reconType == 1  % artifacts chosen by vote
        artifacts = findArtifacts(outfile);
    elseif reconType == 2  % artifacts here are actually normal (included) ones
        artifacts = findArtifacts(outfile);
        artifacts = setdiff(1:size(W), artifacts);
    end
    
    temp_cap = A(:,setdiff(1:size(W,1),artifacts))*St(setdiff(1:size(W,1),artifacts),:);
    [t1,t2,t3] = fileparts(outfile);
    if reconType == 1
        Hdr = writeBESAsb_data(temp_cap, strcat('recon_',t2,'.dat'));
        writeBESAgeneric(strcat('recon_',t2,'.generic'), Hdr);
        writeBESAevt(strcat('recon_',t2,'.evt'), toi, trig);
    elseif reconType == 2
        Hdr = writeBESAsb_data(temp_cap, strcat('reconInverse_',t2,'.dat'));
        writeBESAgeneric(strcat('reconInverse_',t2,'.generic'), Hdr);
        writeBESAevt(strcat('reconInverse_',t2,'.evt'), toi, trig);        
    end

%     blockSize = 800;
%     Hdr.nSamples = 0;
%     Hdr.nChans = size(W,1);
%     Hdr.Fs = 1000;
%     
%     for t = 1:1:size(toi,1)
%        % sources are already demeaned
%        st = St(:, (t-1)*blockSize+1:t*blockSize); 
%        temp_cap = A(:,setdiff(1:size(W,1),artifacts))*st(setdiff(1:size(W,1),artifacts),:);
%        hdr = writeBESAsb_data(temp_cap, strrep(strcat('recon_',outfile),'.out','.dat'));
%        Hdr.nSamples = Hdr.nSamples + hdr.nSamples;
%     end


end
function artifacts = findArtifacts(outfile)
    data = importdata(outfile);
    artifacts = [];
    if isfield(data,'textdata')
        for i = 1:1:size(data.textdata,1)
            ln = strsplit(' ', data.textdata{i});
            if strcmp(ln(2), 'emg') == 1 || strcmp(ln(2), 'peaks') == 1 || strcmp(ln(2), 'eog') == 1 || strcmp(ln(2), 'Reog') == 1
                artifacts = [artifacts str2num(cell2mat(ln(1)))];
            end
        end
    else
        for i = 1:1:size(data,1)
            ln = strsplit(' ', data{i});
            if strcmp(ln(2), 'emg') == 1 || strcmp(ln(2), 'peaks') == 1 || strcmp(ln(2), 'eog') == 1 || strcmp(ln(2), 'Reog') == 1
                artifacts = [artifacts str2num(cell2mat(ln(1)))];
            end
        end
    end
end
function hdr = writeBESAsb_data(X_cap, new_file)
    fid = fopen(new_file, 'w','ieee-le');
    if fid == -1
       fprintf(1,'In writeBESAsb_data:: Can''t WRITE in the folder. Please check your permissions.\n');
       fprintf(1, 'I will quit now... press any key to continue...(no you can''t stop me :D\n');
       pause; 
       exit;
    end
    fwrite(fid, X_cap, 'float32');
    [hdr.nChans,hdr.nSamples] = size(X_cap);
    hdr.Fs = 1000;
    fclose(fid);
end
function writeBESAevt(evtfile, toi, trig)
    % creating event file
    fid = fopen(evtfile,'w');
    if fid == -1
        fprintf(1,'Error creating Event File for Reconstruction of EEG data.\n');
        return;
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    blockSize = 800;
    for t = 1:1:size(toi,1)
       % Following command was used to print trigger at the beginning of
       % each trial.
       % fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*blockSize+1)*1e3, 1, trig(t), strcat('Trial:',num2str(t)));
       % Following command will be used to print trigger at 200 msec after
       % the epoch starts, i.e. where the item was presented.
       fprintf(fid, '%d\t%d\t%d\t%s\n', ((t-1)*blockSize+1+200)*1e3, 1, trig(t), strcat('Trial:',num2str(t)));
    end
    fclose(fid);  
end
function writeBESAgeneric(new_file, hdr)
    fid = fopen(new_file, 'w');
    fprintf(fid,'BESA Generic Data\n\n');
    fprintf(fid,'nChannels=%i\n\n', hdr.nChans);
    fprintf(fid,'sRate=%f\n\n', hdr.Fs);
    fprintf(fid,'nSamples=%i\n\n', hdr.nSamples);
    fprintf(fid,'format=float\n\n');
    fprintf(fid,'file=%s\n\n', strrep(new_file,'.generic','.dat'));
    fclose(fid);
end