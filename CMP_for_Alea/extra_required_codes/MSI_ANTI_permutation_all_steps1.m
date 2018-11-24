function [trialdata cfg timelock]= MSI_ANTI_permutation_all_steps1(datasetname, traillength,evtfile, sfpfile,Fs,prestim,cond,res_trig,elecnum, eeglab_evtfile,EEG_datasetname)

cfg=[];
cfg.dataset =datasetname;
cfg.trialdef.eventtype  = '?';
ft_definetrial(cfg); 
cfg.trialdef.triallength = traillength/Fs; % suppose the length of each trial is 1 sec in MSI and  1.3 in ANTI
cfg = ft_definetrial(cfg);
 
%======Reading EVT file Created by BESA after SOBI

  
 evts = textread(evtfile,'%s');
    tmu = []; 
    trig = [];
    for i = 4:1:size(evts,1)
        if mod(i,3) == 1 
            tmu(end+1) = str2num(cell2mat(evts(i)));
        elseif mod(i,3) == 0
            trig(end+1) = str2num(cell2mat(evts(i)));
        end
    end 

tmu=round(tmu./Fs);  % suppose the Fs=1000Hz and data in evt file is coming in Tmu code trigger order and tmu is in usec
%============================================

%==========Assigning trial (trl) , trial definition and triggers to cfg
trl=[]; cfg.event=[]; 
 % enter the trigger value that you are interested in 
%cond=32;
counter=0;
for i=1:size(tmu,2)
 if ((trig(i)==cond) && (trig(i+1)==res_trig)) 
    counter=counter+1;
    cfg.event(2*counter-1).type='trial';cfg.event(2*counter-1).sample=tmu(i)-prestim+1;cfg.event(2*counter-1).value='[]';
    cfg.event(2*counter).type='backpanel trigger';cfg.event(2*counter).sample=tmu(i)+1;cfg.event(2*counter).value=trig(i);
    trl=[trl;tmu(i)-prestim+1 tmu(i)-prestim+1+(traillength) 0]; % supposed that the epoch data is (-200~+800msec)
 end 
end
cfg.trl=trl;
%====================================================================

% KEEP IN MIND
%I changed the code 'ft_read_data' :
% test whether the requested data segment does not extend over a discontinuous trial boundary
%Iman added
 % hdr.nSamples=1000; hdr.nTrials=1;
%if checkboundary && hdr.nTrials>1  
%==================================================================== 
 
% Pre-procesing  
cfg.trialdef.ntrials     = size(trl,1);  
trialdata = ft_preprocessing(cfg);


% Electrode and SFP file =============================

trialdata.elec = ft_read_sens(sfpfile); % 'C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\ForFieldTrip_EEGLAB\new\recon_6030_anti_stim_-500+1100_forSOBI_ACCA-export-export.sfp');
trialdata.label(1:elecnum)=trialdata.elec.label(4:elecnum+3); % +3 is for Fiducial electrodes 
a=trialdata.elec;
a.label=a.label(4:elecnum+3,:); 
a.elecpos=a.elecpos(4:elecnum+3,:);
a.chanpos=a.chanpos(4:elecnum+3,:);
trialdata.elec=a;
cfg.layout = ft_prepare_layout(cfg, trialdata);
cfg.layout.pos=[cfg.layout.pos(:,2) -cfg.layout.pos(:,1)]; % should od that otherise the plots would not be correct
cfg.showlabels='yes';
trialdata.event=cfg.event;

cfg = [];
cfg.keeptrials = 'yes';
%======================================================
[timelock] = ft_timelockanalysis(cfg, trialdata);
%rmpath('C:\Program Files (x86)\MATLAB\R2011a\toolbox\eeglab9_0_8_6b\external\fieldtrip-partial\fileio');

%=========================================================
% the data in EDF format has oposite polarity!!!!
a=size(trialdata.trial);a=a(2);
for i=1:a
  trialdata.trial{i}=-trialdata.trial{i};
end

