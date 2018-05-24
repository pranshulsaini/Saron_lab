function [ALLEEG EEG] = fieldtrip2eeglab6(dataf)

% load data file ('dataf') preprocessed with fieldtrip
% and show in eeglab viewer

%[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
 ALLEEG=[];
% if exist(dataf,'file')
%   load(dataf)
% end
data=dataf;
% load chanlocs.mat
% EEG.chanlocs = chanlocs; 
EEG.chanlocs = [];

for i=1:size(data.trial,2)
  EEG.data(:,:,i) = single(data.trial{i});
end

% j=0;
% for i=1:size(data.event,2)
%     if  strcmp(data.event(i).type,'backpanel trigger')
%        j=j+1;
%         %event_code=[event_code; data.event(i).sample, data.event(i).type, data.event(i).value];
%        EEG.event(j).sample=(data.event(i).sample-1)/1000;
%        EEG.event(j).type=data.event(i).type;
%        EEG.event(j).value=data.event(i).value;
%     end
% end


EEG.setname    = dataf; %data.cfg.dataset;
EEG.filename   = '';
EEG.filepath   = '';
EEG.subject    = '';
EEG.group      = '';
EEG.condition  = '';
EEG.session    = [];
EEG.comments   = 'preprocessed with fieldtrip';
EEG.nbchan     = size(data.trial{1},1);
EEG.trials     = size(data.trial,2);
EEG.pnts       = size(data.trial{1},2);
EEG.srate      = data.fsample;
EEG.xmin       = data.time{1}(1);
EEG.xmax       = data.time{1}(end);
EEG.times      = data.time{1};
EEG.ref        = []; %'common';
EEG.event      = [];
EEG.epoch      = [];
EEG.icawinv    = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icaact     = [];
EEG.saved      = 'no';

            EEG.icachansind= [];
            EEG.chanlocs= [];
          EEG.urchanlocs= [];
            EEG.chaninfo= [];
             EEG.urevent= [];
    EEG.eventdescription= {};
             EEG.epochdescription= {};
              EEG.reject= [];
               EEG.stats= [];
            EEG.specdata= [];
          EEG.specicaact= [];
          EEG.splinefile= '';
       EEG.icasplinefile= '';
              EEG.dipfit= [];
             EEGhistory= '';
                 EEG.etc= [];


[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
%eeglab(EEG);
%eeglab redraw
%pop_eegplot( EEG, 1, 1, 1);