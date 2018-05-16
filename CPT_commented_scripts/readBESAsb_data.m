% Last modified before Pranshul: 11th April 2012

function data = readBESAsb_data(filename, hdr, beg_sample, end_sample)
    if isempty(findstr(filename,'.dat'))
      filename = [filename,'.dat'];
    end
    % MANISH:: Bytes to skip per sample, is simply equal to 4xnChannels
    bytes_per_sample = 4*hdr.nChans;
    fid=fopen(filename, 'r', 'ieee-le');  % I think it was not necessary. Skipping could have been done without adding 'ieee-le' as well (P)
    % Skippping
    fseek(fid, bytes_per_sample*(beg_sample-1),-1);   % fseek(fileID, offset, origin) (P)  % offset has to be added in terms of bytes (P)
    data=fread(fid,[hdr.nChans,(end_sample-beg_sample+1)],'float32');  % gives output as 88 x n_sample_points
    fclose(fid);
end
 
