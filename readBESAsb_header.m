%(P) means that the comments have been made by Pranshul
% Last modified before Pranshul: 8th Sept 2014
function hdr = readBESAsb_header(filename)

fid=fopen([filename(1:end-4),'.generic'],'r');  % 'end-4' command removed 'dat'
if fid==-1
  fid=fopen([filename(1:end-4),'.gen'],'r');
end


fscanf(fid,'BESA Generic Data\n');    % reads this sentence and then the points goes to the next line (P)
nChannels = fscanf(fid,'nChannels=%i\n');
sRate = fscanf(fid,'sRate=%f\n');
nSamples = fscanf(fid,'nSamples=%i\n');   % they will be less when exported in 'epoch around triggers' format (P)
%format = fscanf(fid,'format=%s');
%file = fscanf(fid,'\nfile=%s');
%prestimulus = fscanf(fid,'prestimulus=%f\n');  % this information is not in my generic file. To get this, you have to export data in 'epoch around triggers' format. I need to comment this line (P)
%epochs = fscanf(fid,'epochs=%i\n'); % this information is not in my generic file. To get this, you have to export data in 'epoch around triggers' format. I need to comment this line (P)
fclose(fid);

% MANISH:: For some weird reason BESA doesn't provide channel labels in the
% DAT file or attached evt/generic file. Thus if we need channel label
% information, which we do in case of FieldTrip, we need to use a fake sfp
% and get the labels - badchannel labels.

% old_sfp = '/Users/mishu/Documents/ShamathaProject/EEG_data/Controls/old_sfp.sfp';
% [new_sfp, drop_list] = select_good_channels(old_sfp, bad_channels);
% fid = fopen(new_sfp,'r');
% labels = {};
% for i = 1:1:nChannels
%     tmp = strsplit(' ', fgets(fid));
%     labels{i} = tmp{1};          %  I think these are only good channel labels (P)
% end

hdr.Fs = sRate;   % hdr is being treated as a structure (P)
hdr.nChans = nChannels ;
hdr.nSamples = nSamples;
%hdr.nSamplesPre = prestimulus; % Will not have this information (P)
%hdr.nTrials = epochs; % will not have this information (P)
% hdr.label = labels;

end