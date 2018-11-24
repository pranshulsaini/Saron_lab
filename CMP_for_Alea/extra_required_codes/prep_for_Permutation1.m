%         cfg = []; 
%         cfg.channel          = {'all'};
%         cfg.latency          = 'all'; 
%    %     cfg.frequency        = foi; 
%         cfg.method           = 'montecarlo'; 
%         cfg.statistic        = 'depSamplesF'; 
%         cfg.correctm         = 'cluster'; 
%         cfg.clusteralpha     = 0.05; 
%         cfg.clusterstatistic = 'maxsum'; 
%         cfg.minnbchan        = 2;   %2 works best.
%         cfg.tail             = 1; 
%         cfg.clustertail      = 1; 
%         cfg.alpha            = 0.05; 
%         cfg.numrandomization = 5000; 
%         
% 
%         subj = endi-starti+1;
%         subj=2;
%         design = zeros(2,2*subj); 
%         for i = 1:1:subj
%            design(1,i) = i;
%            design(1,subj+i) = i;
%            design(1,2*subj+i) = i;
%         end
%         design(2,1:subj) = 1;
%         design(2,subj+1:2*subj) = 2;
%         design(2,2*subj+1:3*subj) = 3;
%         
%         cfg.design = design;
%         cfg.uvar = 1;
%         cfg.ivar = 2;
%         cfg.layout = 'C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\another test\ForFieldtrip_EEGLAB\recon_6030_anti_stim_-500+1100_forSOBI_ACCA_RefFree27Elec_edited.sfp';
%         cfg.elecfile = 'C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\another test\ForFieldtrip_EEGLAB\recon_6030_anti_stim_-500+1100_forSOBI_ACCA_RefFree27Elec_edited.sfp';
%         cfg.neighbourdist = 0.043; % based on BESA
%         cfg.avgoverfreq = 'yes';
%         stat{win} = freqstatistics(cfg, freq_win_GA1{win},freq_win_GA2{win},freq_win_GA3{win});
        
cfg = [];        
cfg.channel          = {'all'};
cfg.latency          = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 1;
cfg.clustertail      = 1;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;
%prepare_neighbours determines what sensors may form clusters
cfg_neighb=[];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, trl1);

design = zeros(1,size(timelock1.trial,1) + size(timelock2.trial,1))% + size(timelock3.trial,1)+size(timelock4.trial,1));
design = zeros(1,size(timelock1.trial,1) + size(timelock2.trial,1)); %+ size(timelock3.trial,1)+size(timelock4.trial,1));

design(1,1:size(timelock1.trial,1)) = 1;
design(1,(size(timelock1.trial,1)+1):(size(timelock1.trial,1) + size(timelock2.trial,1)))= 2;
%design(1,(size(timelock1.trial,1)+1)+(size(timelock1.trial,1) + size(timelock2.trial,1))+1:size(timelock3.trial,1))= 3;
%design(1,(size(timelock1.trial,1)+1)+(size(timelock1.trial,1) + size(timelock2.trial,1))+1+size(timelock3.trial,1)+1:size(timelock4.trial,1))= 4;

cfg.layout = 'C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\another test\ForFieldtrip_EEGLAB\recon_6030_anti_stim_-500+1100_forSOBI_ACCA_RefFree27Elec_edited.sfp';
cfg.elecfile = 'C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\another test\ForFieldtrip_EEGLAB\recon_6030_anti_stim_-500+1100_forSOBI_ACCA_RefFree27Elec_edited.sfp';

cfg.ivar  = 1;                   % number or list with indices, independent variable(s)
%[stat] = ft_timelockstatistics(cfg, timelock1, timelock2,timelock3,timelock4);
[stat] = ft_timelockstatistics(cfg, timelock1, timelock2);


% cfg.layout = ft_prepare_layout(cfg, trl1);
% cfg.layout.pos=[cfg.layout.pos(:,2) -cfg.layout.pos(:,1)]; % should od that otherise the plots would not be correct

cfg = [];
cfg.keeptrials = 'no';  % now keep onl1y the slubject-wise average, not the single trials
avg1 = ft_timelockanalysis(cfg, trl1);
avg2  = ft_timelockanalysis(cfg, trl2);

% Copy the entire timelockanalysis structure to preserve all
% the information it holds in addition to ERP averages.
raweffect     = avg1;
% Then take the difference of the averages.
raweffect.avg = avg1.avg - avg2.avg; 


% Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals = [stat.posclusters(:).prob];
% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% % and now for the negative clusters...
% neg_cluster_pvals = [stat.negclusters(:).prob];
% neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
% neg = ismember(stat.negclusterslabelmat, neg_signif_clust)


timestep = 0.05;		% timestep between time windows for each subplot (in seconds)
sampling_rate = trl1.fsample;	% Data has a temporal resolution of 300 Hz
sample_count = length(stat.time);
					% number of temporal samples in the statistics object
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples


for k = 1:20;
     subplot(4,5,k);
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
     cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel reaches this significance, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).
   
   % Next, check which channels are significant over the
   % entire time interval of interest.
     pos_int = all(pos(:, m(k):m(k+1)), 2);

     cfg.highlight = 'on';
   % Get the index of each significant channel
     cfg.highlightchannel = find(pos_int);
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
     cfg.layout = ft_prepare_layout(cfg, trl1);
     cfg.layout.pos=[cfg.layout.pos(:,2) cfg.layout.pos(:,1)]; % should od that otherise the plots would not be correct
      cfg.showlabels='yes';

     %cfg.layout = 'C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\another test\ForFieldtrip_EEGLAB\recon_6030_anti_stim_-500+1100_forSOBI_ACCA_RefFree27Elec_edited.sfp';
     %cfg.elecfile = 'C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\another test\ForFieldtrip_EEGLAB\recon_6030_anti_stim_-500+1100_forSOBI_ACCA_RefFree27Elec_edited.sfp';
     ft_topoplotER(cfg, raweffect);   
end





















