% This function is for meausring baselines.

function hb = calcIBI(filename)
    % read hdr
    hdr = readBESAsb_header(filename);
    
    % read data
    data = readBESAsb_data(filename, hdr, 1, hdr.nSamples);
    
    % find peaks
    [v i] = findpeaks(-1*data,'MINPEAKDISTANCE',1200,'sortstr','descend');
    i = sort(i);
    
    % calculate ibi on the whole
    hbibi = diff(i);
    
    % remove artifacts
    
    % find IBI per condition/trial
    
    

end

function hdr = readBESAsb_header(filename)

    fid=fopen([filename(1:end-4),'.generic'],'r');
    if fid==-1
      fid=fopen([filename(1:end-4),'.gen'],'r');
    end

    fscanf(fid,'BESA Generic Data\n');
    nChannels = fscanf(fid,'nChannels=%i\n');
    sRate = fscanf(fid,'sRate=%f\n');
    nSamples = fscanf(fid,'nSamples=%i\n');
    format = fscanf(fid,'format=%s');
    file = fscanf(fid,'\nfile=%s');
    prestimulus = fscanf(fid,'prestimulus=%f\n');
    epochs = fscanf(fid,'epochs=%i\n');
    fclose(fid);

    hdr.Fs = sRate;
    hdr.nChans = nChannels ;
    hdr.nSamples = nSamples;
    hdr.nSamplesPre = prestimulus;
    hdr.nTrials = epochs;

end
function data = readBESAsb_data(filename, hdr, beg_sample, end_sample)
    if isempty(findstr(filename,'.dat'))
      filename = [filename,'.dat'];
    end
    % MANISH:: Bytes to skip per sample, is simply equal to
    % 4xnChannels
    bytes_per_sample = 4*hdr.nChans;
    fid=fopen(filename, 'r', 'ieee-le');
    % Skippping
    fseek(fid, bytes_per_sample*(beg_sample-1),-1);
    data=fread(fid,[hdr.nChans,(end_sample-beg_sample+1)],'float32');
    fclose(fid);
end

