%% Plots GFP as demanded by Cliff
clc;

%base_dir = 'C:\Users\plsaini\Box Sync\Stroop\STR_SOBI_Output\';
base_dir = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Filtered\';

sub_ID = 125;
subject_ID = sub_ID;  % I have used these two terms interchangebly
type = 'retreat';

% local_dir_pre = strcat('STR',num2str(subject_ID),'11_aux_NoBad_AvgRef_forSOBI\STR',num2str(subject_ID),'11_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR',num2str(subject_ID),'11_aux_NoBad_AvgRef_forSOBI_recon\STR');
% 
% local_dir_mid = strcat('STR',num2str(subject_ID),'12_aux_NoBad_AvgRef_forSOBI\STR',num2str(subject_ID),'12_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR',num2str(subject_ID),'12_aux_NoBad_AvgRef_forSOBI_recon\STR');
% 
% local_dir_post = strcat('STR',num2str(subject_ID),'13_aux_NoBad_AvgRef_forSOBI\STR',num2str(subject_ID),'13_aux_NoBad_AvgRef_forSOBI_Sources_STR\STR',num2str(subject_ID),'13_aux_NoBad_AvgRef_forSOBI_recon\STR');

% % data at pre-assessment
data_pre_neu =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'11_neu_corr.Export.Export.eph'));
data_pre_incong =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'11_incong_corr.Export.Export.eph'));
data_pre_cong =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'11_cong_corr.Export.Export.eph'));

% finding max value at pre-assessment
max_pre_neu = max(data_pre_neu);
max_pre_neu = max_pre_neu(1);
max_pre_incong = max(data_pre_incong);
max_pre_incong = max_pre_incong(1);
max_pre_cong = max(data_pre_cong);
max_pre_cong = max_pre_cong(1);
max_pre = max([max_pre_neu,max_pre_incong,max_pre_cong]);


% data at mid-assessment
data_mid_neu =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'12_neu_corr.Export.Export.eph'));
data_mid_incong =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'12_incong_corr.Export.Export.eph'));
data_mid_cong =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'12_cong_corr.Export.Export.eph'));

% finding max value at mid-assessment
max_mid_neu = max(data_mid_neu);
max_mid_neu = max_mid_neu(1);
max_mid_incong = max(data_mid_incong);
max_mid_incong = max_mid_incong(1);
max_mid_cong = max(data_mid_cong);
max_mid_cong = max_mid_cong(1);
max_mid = max([max_mid_neu,max_mid_incong,max_mid_cong]);

% data at post-assessment 
data_post_neu =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'13_neu_corr.Export.Export.eph'));
data_post_incong =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'13_incong_corr.Export.Export.eph'));
data_post_cong =  importdata(strcat(base_dir,'STR',num2str(subject_ID),'13_cong_corr.Export.Export.eph'));

% finding max value at post-assessment
max_post_neu = max(data_post_neu);
max_post_neu = max_post_neu(1);
max_post_incong = max(data_post_incong);
max_post_incong = max_post_incong(1);
max_post_cong = max(data_post_cong);
max_post_cong = max_post_cong(1);
max_post = max([max_post_neu,max_post_incong,max_post_cong]);

super_max = max([max_pre, max_mid, max_post]);
super_max = fix(super_max) +1; % it will be used to fix y-axis length

time = linspace(-100,1200,2663);  % 2663 points between -100 ms and 1200 ms

% plotting at pre-retreat
FigH = figure('Position', get(0, 'Screensize')); % will help in saving picture in a full screen size
subplot(3,1,1)
plot(time, data_pre_neu(308:2970,1),'k');
hold on;
plot(time, data_pre_incong(308:2970,1),'r');
plot(time, data_pre_cong(308:2970,1),'g');
lgd = legend({'neu','incong','cong'},'Location','northwest');
lgd.FontSize = 7;
title(strcat(type,'-',num2str(sub_ID),'-pre'));
ylim([0 super_max]);

% plotting at mid-retreat
subplot(3,1,2)
plot(time, data_mid_neu(308:2970,1),'k');
hold on;
plot(time, data_mid_incong(308:2970,1),'r');
plot(time, data_mid_cong(308:2970,1),'g');
lgd = legend({'neu','incong','cong'},'Location','northwest');
lgd.FontSize = 7;
title(strcat(type,'-',num2str(sub_ID),'-mid'));
ylim([0 super_max]);

% plotting at mid-retreat
subplot(3,1,3)
plot(time, data_post_neu(308:2970,1),'k');
hold on;
plot(time, data_post_incong(308:2970,1),'r');
plot(time, data_post_cong(308:2970,1),'g');
lgd = legend({'neu','incong','cong'},'Location','northwest');
lgd.FontSize = 7;
title(strcat(type,'-',num2str(sub_ID),'-post'));
xlabel ('time (ms)')
ylim([0 super_max]);

saveas(gcf,strcat('C:\Users\plsaini\Box Sync\Stroop\GFP\','STR',num2str(sub_ID),'GFP'),'jpg');
