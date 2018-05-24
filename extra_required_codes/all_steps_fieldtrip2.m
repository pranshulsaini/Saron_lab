cfg=[];
% Enter the path for +EDF file 
%cfg.dataset ='C:\Users\imanmr\Desktop\Fieldtrip_test\complete data2\recon_ASDMSI_021_forSOBI_ACCA-export.edf';
cfg.dataset ='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_021\Goto SOBI\SOBI_Output\ASDMSI_021_forSOBI\forFieldtrip_EEGLAB\recon_ASDMSI_021_forSOBI_ACCA-export.edf';
cfg.trialdef.eventtype  = '?';
ft_definetrial(cfg); 

cfg.trialdef.triallength = 1; % suppose the length of each trial is 1 sec
cfg = ft_definetrial(cfg);
 
%======Reading EVT file Created by BESA after SOBI
% evt file should have event time, code and trigger number 
%evtfile='C:\Users\imanmr\Desktop\Fieldtrip_test\complete data2\recon_ASDMSI_021_forSOBI_ACCA-export.evt';
evtfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_021\Goto SOBI\SOBI_Output\ASDMSI_021_forSOBI\forFieldtrip_EEGLAB\recon_ASDMSI_021_forSOBI_ACCA.evt';

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

tmu=round(tmu./1000);  % suppoese the Fs=1000Hz and data in evt file is coming in Tmu code trigger order and tmu is in usec
%============================================

%==========Assigning trial (trl) , trial definition and triggers to cfg
trl=[]; cfg.event=[]; 
 % enter the condition value that you are interested in 
cond=16;
counter=0;
for i=1:size(tmu,2)
 if trig(i)==cond
    counter=counter+1;
    cfg.event(2*counter-1).type='trial';cfg.event(2*counter-1).sample=tmu(i)-199;cfg.event(2*counter-1).value='[]';
    cfg.event(2*counter).type='backpanel trigger';cfg.event(2*counter).sample=tmu(i)+1;cfg.event(2*counter).value=trig(i);
    trl=[trl;tmu(i)-199 tmu(i)-199+1000 0]; % supposed that the epoch data is (-200~+800msec)
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
%trialdata.elec = ft_read_sens('C:\Users\imanmr\Desktop\Fieldtrip_test\average\average_tst._2cond.sfp');
% Enter sfp file here
%trialdata.elec = ft_read_sens('C:\Users\imanmr\Desktop\Fieldtrip_test\complete data2\recon_ASDMSI_021_forSOBI_ACCA-export.sfp');
trialdata.elec = ft_read_sens('C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_021\Goto SOBI\SOBI_Output\ASDMSI_021_forSOBI\forFieldtrip_EEGLAB\recon_ASDMSI_021_forSOBI_ACCA-export.sfp');
trialdata.label(1:132)=trialdata.elec.label(4:135);
a=trialdata.elec;
a.label=a.label(4:135,:);
a.elecpos=a.elecpos(4:135,:);
a.chanpos=a.chanpos(4:135,:);
trialdata.elec=a;
cfg.layout = ft_prepare_layout(cfg, trialdata);
cfg.layout.pos=[cfg.layout.pos(:,2) -cfg.layout.pos(:,1)]; % should od that otherise the plots would not be correct
cfg.showlabels='yes';
trialdata.event=cfg.event;
%===================================================
[timelock] = ft_timelockanalysis(cfg, trialdata);
%timelock.avg=-timelock.avg; % data looks inverted by comparing to BESA
cfg2 = ft_databrowser(cfg,timelock); 
cfg = ft_databrowser(cfg,trialdata);   