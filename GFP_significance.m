%% GFP plots for significant time windows. The significant time windows are extracted from Ragu

% Also we need grand average GFP waveform overlays of the two groups by condition separately by assessment, with significant time periods of difference. 
% Possibly also overlays of 6 waveforms - by condition and group at each assessment.

results = importdata('C:\Users\plsaini\Box Sync\Stroop\Ragu\GFP.txt');

p_cong = results.data(:,2);
p_assess = results.data(:,3);
p_cong_assess = results.data(:,4);
p_group = results.data(:,5);
p_cong_group = results.data(:,6);
p_assess_group = results.data(:,7);
p_cong_assess_group = results.data(:,8);

% this will store the points which have p_values for more than 20 ms (41 points)
index_1 = [];
index_2 = [];
index_3 = [];
index_4 = [];
index_5 = [];
index_6 = [];
index_7 = [];

track = ones(8,1);
for i = 1:2971-40
    score =  zeros(8,1);
    for j = 0:40
        if (p_cong(i+j)<0.05)
            score(1) = score(1) + 1;
        end
        
        if (p_assess(i+j)<0.05)
            score(2) = score(2) + 1;
        end
        
        if (p_cong_assess(i+j)<0.05)
            score(3) = score(3) + 1;
        end
        
        if (p_group(i+j)<0.05)
            score(4) = score(4) + 1;
        end
        
        if (p_cong_group(i+j)<0.05)
            score(5) = score(5) + 1;
        end
    
        if (p_assess_group(i+j)<0.05)
            score(6) = score(6) + 1;
        end

        if (p_cong_assess_group(i+j)<0.05)
            score(7) = score(7) + 1;
        end
    end
    
    if (score(1)==41)
        index_1(track(1),1) = i;
        index_1(track(1),2) = i+40;
        track(1) = track(1) + 1;
    end
    
    if (score(2)==41)
        index_2(track(2),1) = i;
        index_2(track(2),2) = i+40;
        track(2) = track(2) + 1;
    end
    
    if (score(3)==41)
        index_3(track(3),1) = i;
        index_3(track(3),2) = i+40;
        track(3) = track(3) + 1;
    end
    
    if (score(4)==41)
        index_4(track(4),1) = i;
        index_4(track(4),2) = i+40;
        track(4) = track(4) + 1;
    end
    
    if (score(5)==41)
        index_5(track(5),1) = i;
        index_5(track(5),2) = i+40;
        track(5) = track(5) + 1;
    end
    
    if (score(6)==41)
        index_6(track(6),1) = i;
        index_6(track(6),2) = i+40;
        track(6) = track(6) + 1;
    end
    
    if (score(7)==41)
        index_7(track(7),1) = i;
        index_7(track(7),2) = i+40;
        track(7) = track(7) + 1;
    end
    
end


%significant intervals
sig_time = [281:333, 935:985, 1008:1067, 1162:1364, 1371:1617, 2068:2148, 2167:2688,  2802:2970 ];
sig_time_ms = sig_time*0.48828125 - 250;

%% Overlay of two groups by condition. Plot: 3 x 1
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Neu\Avg STR All\';
neu_control = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
neu_control_sig = neu_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Cong\Avg STR All\';
cong_control = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
cong_control_sig = cong_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Incong\Avg STR All\';
incong_control = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
incong_control_sig = incong_control(sig_time);

file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Neu\Avg STR All\';
neu_retreat = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
neu_retreat_sig = neu_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Cong\Avg STR All\';
cong_retreat = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
cong_retreat_sig = cong_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Incong\Avg STR All\';
incong_retreat = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
incong_retreat_sig = incong_retreat(sig_time);

super_max = max([max(neu_control), max(cong_control), max(incong_control), max(neu_retreat), max(cong_retreat), max(incong_retreat)]);

time = linspace(-250,1200,2970);  % 2970 points between -250 ms and 1200 ms


FigH = figure('Position', get(0, 'Screensize'));%,'visible','off'); % will help in saving picture in a full screen size. Visible-off feature keeps figures from popping up
subplot(3,1,1);
plot(time, neu_control,'r');
hold on;
plot(time, neu_retreat,'b');
title(strcat('Neutral'));
ylim([0 super_max]);
plot(sig_time_ms, neu_control_sig,'r.');
plot(sig_time_ms, neu_retreat_sig,'b.');
lgd = legend({'control','retreat','controlSig','retreatSig'},'Location','northwest');
lgd.FontSize = 7;

subplot(3,1,2);
plot(time, cong_control,'r');
hold on;
plot(time, cong_retreat,'b');
title(strcat('Congruent'));
ylim([0 super_max]);
plot(sig_time_ms, cong_control_sig,'r.');
plot(sig_time_ms, cong_retreat_sig,'b.');
lgd = legend({'control','retreat','controlSig','retreatSig'},'Location','northwest');
lgd.FontSize = 7;

subplot(3,1,3);
plot(time, incong_control,'r');
hold on;
plot(time, incong_retreat,'b');
title(strcat('Incongruent'));
ylim([0 super_max]);
plot(sig_time_ms, incong_control_sig,'r.');
plot(sig_time_ms, incong_retreat_sig,'b.');
lgd = legend({'control','retreat','controlSig','retreatSig'},'Location','northwest');
lgd.FontSize = 7;

saveas(gcf,strcat('C:\Users\plsaini\Box Sync\Stroop\GFP\','GroupxCondition'),'jpg');


%% Overlay of two groups by assessment. Plot: 3 x 1

file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Pre\Avg STR All\';
pre_control = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
pre_control_sig = pre_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Mid\Avg STR All\';
mid_control = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
mid_control_sig = mid_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Post\Avg STR All\';
post_control = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
post_control_sig = post_control(sig_time);

file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Pre\Avg STR All\';
pre_retreat = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
pre_retreat_sig = pre_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Mid\Avg STR All\';
mid_retreat = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
mid_retreat_sig = mid_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Post\Avg STR All\';
post_retreat = importdata(strcat(file_loc,'Avg STR All.Export.ep'));
post_retreat_sig = post_retreat(sig_time);

super_max = max([max(pre_control), max(mid_control), max(post_control), max(pre_retreat), max(mid_retreat), max(post_retreat)]);

time = linspace(-250,1200,2970);  % 2970 points between -250 ms and 1200 ms


FigH = figure('Position', get(0, 'Screensize'));%,'visible','off'); % will help in saving picture in a full screen size. Visible-off feature keeps figures from popping up
subplot(3,1,1);
plot(time, pre_control,'r');
hold on;
plot(time, pre_retreat,'b');
title(strcat('Pre'));
ylim([0 super_max]);
plot(sig_time_ms, pre_control_sig,'r.');
plot(sig_time_ms, pre_retreat_sig,'b.');
lgd = legend({'control','retreat','controlSig','retreatSig'},'Location','northwest');
lgd.FontSize = 7;

subplot(3,1,2);
plot(time, mid_control,'r');
hold on;
plot(time, mid_retreat,'b');
title(strcat('Mid'));
ylim([0 super_max]);
plot(sig_time_ms, mid_control_sig,'r.');
plot(sig_time_ms, mid_retreat_sig,'b.');
lgd = legend({'control','retreat','controlSig','retreatSig'},'Location','northwest');
lgd.FontSize = 7;

subplot(3,1,3);
plot(time, post_control,'r');
hold on;
plot(time, post_retreat,'b');
title(strcat('Post'));
ylim([0 super_max]);
plot(sig_time_ms, post_control_sig,'r.');
plot(sig_time_ms, post_retreat_sig,'b.');
lgd = legend({'control','retreat','controlSig','retreatSig'},'Location','northwest');
lgd.FontSize = 7;

saveas(gcf,strcat('C:\Users\plsaini\Box Sync\Stroop\GFP\','GroupxAssessment'),'jpg');


%% overlays of 6 waveforms - by condition and group at each assessment

%pre_control
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Pre_Neu\Avg StrControl All\';
pre_neu_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
pre_neu_control_sig = pre_neu_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Pre_Cong\Avg StrControl All\';
pre_cong_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
pre_cong_control_sig = pre_cong_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Pre_Incong\Avg StrControl All\';
pre_incong_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
pre_incong_control_sig = pre_incong_control(sig_time);

%mid_control
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Mid_Neu\Avg StrControl All\';
mid_neu_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
mid_neu_control_sig = mid_neu_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Mid_Cong\Avg StrControl All\';
mid_cong_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
mid_cong_control_sig = mid_cong_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Mid_Incong\Avg StrControl All\';
mid_incong_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
mid_incong_control_sig = mid_incong_control(sig_time);

%post_control
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Post_Neu\Avg StrControl All\';
post_neu_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
post_neu_control_sig = post_neu_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Post_Cong\Avg StrControl All\';
post_cong_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
post_cong_control_sig = post_cong_control(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Control\Post_Incong\Avg StrControl All\';
post_incong_control = importdata(strcat(file_loc,'Avg StrControl All.Export.ep'));
post_incong_control_sig = post_incong_control(sig_time);


%pre_retreat
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Pre_Neu\Avg StrRetreat All\';
pre_neu_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
pre_neu_retreat_sig = pre_neu_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Pre_Cong\Avg StrRetreat All\';
pre_cong_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
pre_cong_retreat_sig = pre_cong_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Pre_Incong\Avg StrRetreat All\';
pre_incong_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
pre_incong_retreat_sig = pre_incong_retreat(sig_time);

%mid_retreat
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Mid_Neu\Avg StrRetreat All\';
mid_neu_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
mid_neu_retreat_sig = mid_neu_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Mid_Cong\Avg StrRetreat All\';
mid_cong_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
mid_cong_retreat_sig = mid_cong_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Mid_Incong\Avg StrRetreat All\';
mid_incong_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
mid_incong_retreat_sig = mid_incong_retreat(sig_time);

%post_retreat
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Post_Neu\Avg StrRetreat All\';
post_neu_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
post_neu_retreat_sig = post_neu_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Post_Cong\Avg StrRetreat All\';
post_cong_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
post_cong_retreat_sig = post_cong_retreat(sig_time);
file_loc = 'C:\Users\plsaini\Box Sync\Stroop\Cartool\Grand_Avg\Retreat\Post_Incong\Avg StrRetreat All\';
post_incong_retreat = importdata(strcat(file_loc,'Avg StrRetreat All.Export.ep'));
post_incong_retreat_sig = post_incong_retreat(sig_time);


super_max = max([max(pre_neu_control),max(mid_neu_control),max(post_neu_control), max(pre_cong_control),max(mid_cong_control), max(post_cong_control), max(pre_incong_control),max(mid_incong_control), max(post_incong_control), max(pre_neu_retreat),max(mid_neu_retreat),max(post_neu_retreat), max(pre_cong_retreat),max(mid_cong_retreat), max(post_cong_retreat), max(pre_incong_retreat),max(mid_incong_retreat), max(post_incong_retreat)]);

time = linspace(-250,1200,2970);  % 2970 points between -250 ms and 1200 ms


FigH = figure('Position', get(0, 'Screensize'));%,'visible','off'); % will help in saving picture in a full screen size. Visible-off feature keeps figures from popping up
subplot(3,1,1);
h1 = plot(time, pre_neu_control,'b');
hold on;
h2 = plot(time, pre_neu_retreat,'k');
h3 = plot(time, pre_cong_control,'r');
h4 = plot(time, pre_cong_retreat,'g');
h5 = plot(time, pre_incong_control,'c');
h6 = plot(time, pre_incong_retreat,'m');

title(strcat('Pre'));
ylim([0 super_max]);
plot(sig_time_ms, pre_neu_control_sig,'b.');
plot(sig_time_ms, pre_neu_retreat_sig,'k.');
plot(sig_time_ms, pre_cong_control_sig,'r.');
plot(sig_time_ms, pre_cong_retreat_sig,'g.');
plot(sig_time_ms, pre_incong_control_sig,'c.');
plot(sig_time_ms, pre_incong_retreat_sig,'m.');
lgd = legend([h1 h2 h3 h4 h5 h6],{'neuControl','neuRetreat','congControl','congRetreat','incongControl','incongRetreat'},'Location','northwest');
lgd.FontSize = 7;

subplot(3,1,2);
h1 = plot(time, mid_neu_control,'b');
hold on;
h2 = plot(time, mid_neu_retreat,'k');
h3 = plot(time, mid_cong_control,'r');
h4 = plot(time, mid_cong_retreat,'g');
h5 = plot(time, mid_incong_control,'c');
h6 = plot(time, mid_incong_retreat,'m');

title(strcat('Mid'));
ylim([0 super_max]);
plot(sig_time_ms, mid_neu_control_sig,'b.');
plot(sig_time_ms, mid_neu_retreat_sig,'k.');
plot(sig_time_ms, mid_cong_control_sig,'r.');
plot(sig_time_ms, mid_cong_retreat_sig,'g.');
plot(sig_time_ms, mid_incong_control_sig,'c.');
plot(sig_time_ms, mid_incong_retreat_sig,'m.');
lgd = legend([h1 h2 h3 h4 h5 h6],{'neuControl','neuRetreat','congControl','congCetreat','incongControl','incongRetreat','retreatSig'},'Location','northwest');
lgd.FontSize = 7;

subplot(3,1,3);
h1 = plot(time, post_neu_control,'b');
hold on;
h2 = plot(time, post_neu_retreat,'k');
h3 = plot(time, post_cong_control,'r');
h4 = plot(time, post_cong_retreat,'g');
h5 = plot(time, post_incong_control,'c');
h6 = plot(time, post_incong_retreat,'m');

title(strcat('Post'));
ylim([0 super_max]);
plot(sig_time_ms, post_neu_control_sig,'b.');
plot(sig_time_ms, post_neu_retreat_sig,'k.');
plot(sig_time_ms, post_cong_control_sig,'r.');
plot(sig_time_ms, post_cong_retreat_sig,'g.');
plot(sig_time_ms, post_incong_control_sig,'c.');
plot(sig_time_ms, post_incong_retreat_sig,'m.');
lgd = legend([h1 h2 h3 h4 h5 h6],{'neuControl','neuRetreat','congControl','congRetreat','incongControl','incongRetreat'},'Location','northwest');
lgd.FontSize = 7;

saveas(gcf,strcat('C:\Users\plsaini\Box Sync\Stroop\GFP\','GroupxAssessmentxCondition'),'jpg');
