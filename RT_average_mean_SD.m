%% This file will average reaction time data across subjects for each group, condition, and assessment. Keep in mind that all 128 trials have been used

clc;
clear all;
%online link: https://www.mathworks.com/matlabcentral/answers/306876-how-do-i-read-only-a-specific-line-while-reading-a-text-file-in-matlab
%online link: https://www.mathworks.com/matlabcentral/answers/339342-how-to-read-a-letter-from-a-string-in-a-cell

group = 'Retreat';

sub_IDs = importdata('C:\Users\plsaini\Box Sync\Stroop\Subject_IDs.txt');
num_sub = size(sub_IDs,1);
cond = zeros(num_sub,1);

% 8 condition sequences. There are eight different orders of sequences and each has been assigned four columns in the excel file. The first column is
% the condition{-1,0,-1} for trials, the second one is the block number [1-16], the third one as block type {-1,0,1}, and the fourth one as blank col
conditions = importdata('C:\Users\plsaini\Box Sync\Stroop\conditions_sorted.xlsx');
conditions = cell2mat(struct2cell(conditions));

%arrays to store avg RT
RT_avg_pre = zeros(round(num_sub/3),5); % 5 conditions
RT_avg_mid = zeros(round(num_sub/3),5); % 5 conditions
RT_avg_post = zeros(round(num_sub/3),5); % 5 conditions

%arrays to store RT std 
RT_sd_pre = zeros(round(num_sub/3),5); % 5 conditions
RT_sd_mid = zeros(round(num_sub/3),5); % 5 conditions
RT_sd_post = zeros(round(num_sub/3),5); % 5 conditions

count_pre = 0;
count_mid = 0;
count_post = 0;


for i = 1:num_sub % 2nd file has incomplete event file. So a mismatch arises b/w log and event file
    i
    counter  = 0;
    
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\Logfiles\',sub_IDs{i}); % log file
    fid=fopen(filename);
    linenum = 7;      % This line number has condition written
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    cond_num = C{1}{1}(29);  % finding condition number
    cond(i) = str2num(cond_num);
    fclose(fid);
    
    
    %log_file
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\Logfiles\',sub_IDs{i}); 
    fid = fopen(filename);
    log_file = textscan(fid,'%d %s %s %s %d %d','headerlines',98);
    stim = log_file{1,3}(:);
    stim_len = cellfun('length',stim);
    resp = log_file{1,4}(:);
    resp_len = cellfun('length', resp);
    outcome = stim_len == resp_len;
    outcome = outcome + 399; % 400 for correct, 399 for incorrect
    fclose(fid);
    
    % creating different markers for different conditions
    cond_seq = (conditions(:,4*cond(i)-3) + conditions(:,4*cond(i)-1) + 3)*1000; % the sum of 1st two terms lead to -2(incongruent), -1 (neutral in incongruent), 0 (neutral in neutral)
    %, 1 (neutral in congruent), 2 congruent. After adding 3, now it will be 1000,2000,3000,4000,5000 
    
    % RT averageing is being done for pre and only for correct trials
    if (sub_IDs{i}(8) == '1')
        count_pre = count_pre + 1;
        RT_incong = [];
        RT_incong_neu = [];
        RT_neu = [];
        RT_cong_neu = [];
        RT_cong = [];
        for j = 1:256
            if ((cond_seq(j) == 1000) && (outcome(j) == 400))
                RT_incong(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 2000) && (outcome(j) == 400))
                RT_incong_neu(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 3000) && (outcome(j) == 400))
                RT_neu(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 4000) && (outcome(j) == 400))
                RT_cong_neu(end + 1) = log_file{1,5}(j);    
            elseif ((cond_seq(j) == 5000) && (outcome(j) == 400))
                RT_cong(end + 1) = log_file{1,5}(j);
            elseif (outcome(j) == 400)   % we will be checking only for correct trials
                display('NO Condition matched')
            end            
        end
    
        
    RT_avg_pre(count_pre,:) =[mean(RT_incong), mean(RT_incong_neu), mean(RT_neu), mean(RT_cong_neu), mean(RT_cong)];
    RT_sd_pre(count_pre,:) =[std(RT_incong), std(RT_incong_neu), std(RT_neu), std(RT_cong_neu), std(RT_cong)];
    counter = 1;
    end
    
    
    % RT averageing is being done for mid and only for correct trials
    if (sub_IDs{i}(8) == '2')
        count_mid = count_mid + 1;
        RT_incong = [];
        RT_incong_neu = [];
        RT_neu = [];
        RT_cong_neu = [];
        RT_cong = [];
        for j = 1:256
            if ((cond_seq(j) == 1000) && (outcome(j) == 400))
                RT_incong(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 2000) && (outcome(j) == 400))
                RT_incong_neu(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 3000) && (outcome(j) == 400))
                RT_neu(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 4000) && (outcome(j) == 400))
                RT_cong_neu(end + 1) = log_file{1,5}(j);    
            elseif ((cond_seq(j) == 5000) && (outcome(j) == 400))
                RT_cong(end + 1) = log_file{1,5}(j);
            elseif (outcome(j) == 400)   % we will be checking only for correct trials
                display('NO Condition matched')
            end            
        end
        
    RT_avg_mid(count_mid,:) =[mean(RT_incong), mean(RT_incong_neu), mean(RT_neu), mean(RT_cong_neu), mean(RT_cong)];
    RT_sd_mid(count_mid,:) =[std(RT_incong), std(RT_incong_neu), std(RT_neu), std(RT_cong_neu), std(RT_cong)];
    counter = 1;
    end
    
    
    % RT averageing is being done for post and only for correct trials
    if (sub_IDs{i}(8) == '3')
        count_post = count_post + 1;
        RT_incong = [];
        RT_incong_neu = [];
        RT_neu = [];
        RT_cong_neu = [];
        RT_cong = [];
        for j = 1:256
            if ((cond_seq(j) == 1000) && (outcome(j) == 400))
                RT_incong(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 2000) && (outcome(j) == 400))
                RT_incong_neu(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 3000) && (outcome(j) == 400))
                RT_neu(end + 1) = log_file{1,5}(j);
            elseif ((cond_seq(j) == 4000) && (outcome(j) == 400))
                RT_cong_neu(end + 1) = log_file{1,5}(j);    
            elseif ((cond_seq(j) == 5000) && (outcome(j) == 400))
                RT_cong(end + 1) = log_file{1,5}(j);
            elseif (outcome(j) == 400)   % we will be checking only for correct trials
                display('NO Condition matched')
            end            
        end
        
    RT_avg_post(count_post,:) =[mean(RT_incong), mean(RT_incong_neu), mean(RT_neu), mean(RT_cong_neu), mean(RT_cong)];
    RT_sd_post(count_post,:) =[std(RT_incong), std(RT_incong_neu), std(RT_neu), std(RT_cong_neu), std(RT_cong)];
    counter = 1;
    end
    
    if(counter == 0)
        display('Not pre, not mid, not post!!!!!! CHECK!!')
    end
end

super_max = max([max(RT_avg_pre(:)),max(RT_avg_mid(:)),max(RT_avg_post(:))]);
super_max = max(super_max);

super_min = min([min(RT_avg_pre(:)),min(RT_avg_mid(:)),min(RT_avg_post(:))]);
super_min = min(super_min);

%%%%%%%%%%%%%%%%%%%%%%saving data%%%%%%%%%%%%%%%%%%%

save(strcat(group,'_RT_avg_pre_mean'),'RT_avg_pre');
save(strcat(group,'_RT_avg_mid_mean'),'RT_avg_mid');
save(strcat(group,'_RT_avg_post_mean'),'RT_avg_post');

save(strcat(group,'_RT_sd_pre_mean'),'RT_sd_pre');
save(strcat(group,'_RT_sd_mid_mean'),'RT_sd_mid');
save(strcat(group,'_RT_sd_post_mean'),'RT_sd_post');

% averaging across subjects as well
Grand_RT_avg_pre = mean(RT_avg_pre);
Grand_RT_avg_mid = mean(RT_avg_mid);
Grand_RT_avg_post = mean(RT_avg_post);

save(strcat(group,'_Grand_RT_avg_pre_mean'),'Grand_RT_avg_pre');
save(strcat(group,'_Grand_RT_avg_mid_mean'),'Grand_RT_avg_mid');
save(strcat(group,'_Grand_RT_avg_post_mean'),'Grand_RT_avg_post');

% averaging across subjects as well
Grand_RT_sd_pre = mean(RT_sd_pre);
Grand_RT_sd_mid = mean(RT_sd_mid);
Grand_RT_sd_post = mean(RT_sd_post);

save(strcat(group,'_Grand_RT_sd_pre_mean'),'Grand_RT_sd_pre');
save(strcat(group,'_Grand_RT_sd_mid_mean'),'Grand_RT_sd_mid');
save(strcat(group,'_Grand_RT_sd_post_mean'),'Grand_RT_sd_post');


%% assessment wise plotting
figure;
subplot(3,1,1);
RT_avg_pre3 = RT_avg_pre(:,1:2:end); 
plot(RT_avg_pre3);
%legend({'incong','incong-neu','neu','cong-neu','cong'},'Location','southeast')
legend({'incong','neu','cong'},'Location','southeast')
title(strcat(group,' RT at pre'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

subplot(3,1,2);
RT_avg_mid3 = RT_avg_mid(:,1:2:end);
plot(RT_avg_mid3);
%legend({'incong','incong-neu','neu','cong-neu','cong'},'Location','southeast')
legend({'incong','neu','cong'},'Location','southeast')
title(strcat(group,' RT at mid'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

subplot(3,1,3)
RT_avg_post3 = RT_avg_post(:,1:2:end);
plot(RT_avg_post3);
%legend({'incong','incong-neu','neu','cong-neu','cong'},'Location','southeast')
legend({'incong','neu','cong'},'Location','southeast')
title(strcat(group,' RT at post'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

%% condition-wise plotting
figure;
subplot(3,1,1);
plot(RT_avg_pre(:,1));
hold on;
plot(RT_avg_mid(:,1));
plot(RT_avg_post(:,1));
legend({'pre','mid','post'},'Location','southeast')
title(strcat(group,' RT: incongruent'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

subplot(3,1,2);
plot(RT_avg_pre(:,3));
hold on;
plot(RT_avg_mid(:,3));
plot(RT_avg_post(:,3));
legend({'pre','mid','post'},'Location','southeast')
title(strcat(group,' RT: neutral'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);


subplot(3,1,3);
plot(RT_avg_pre(:,5));
hold on;
plot(RT_avg_mid(:,5));
plot(RT_avg_post(:,5));
legend({'pre','mid','post'},'Location','southeast')
title(strcat(group,' RT: congruent'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);




%% plotting grand RT averages with standard deviation
Retreat_pre_RT = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_avg_pre_mean.mat');
Retreat_mid_RT = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_avg_mid_mean.mat');
Retreat_post_RT = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_avg_post_mean.mat');
Control_pre_RT = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_avg_pre_mean.mat');
Control_mid_RT = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_avg_mid_mean.mat');
Control_post_RT = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_avg_post_mean.mat');

data_Retreat_pre = Retreat_pre_RT.Grand_RT_avg_pre;
data_Retreat_mid = Retreat_mid_RT.Grand_RT_avg_mid;
data_Retreat_post = Retreat_post_RT.Grand_RT_avg_post;
data_Control_pre = Control_pre_RT.Grand_RT_avg_pre;
data_Control_mid = Control_mid_RT.Grand_RT_avg_mid;
data_Control_post = Control_post_RT.Grand_RT_avg_post;

data_Control = [data_Control_pre', data_Control_mid', data_Control_post'];
data_Retreat = [data_Retreat_pre', data_Retreat_mid', data_Retreat_post'];

super_max = max([max(data_Control(:)), max(data_Retreat(:))]);
super_min = min([min(data_Control(:)), min(data_Retreat(:))]);

Retreat_pre_sd = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_sd_pre_mean.mat');
Retreat_mid_sd = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_sd_mid_mean.mat');
Retreat_post_sd = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_sd_post_mean.mat');
Control_pre_sd = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_sd_pre_mean.mat');
Control_mid_sd = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_sd_mid_mean.mat');
Control_post_sd = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_sd_post_mean.mat');

sd_Retreat_pre = Retreat_pre_sd.Grand_RT_sd_pre;
sd_Retreat_mid = Retreat_mid_sd.Grand_RT_sd_mid;
sd_Retreat_post = Retreat_post_sd.Grand_RT_sd_post;
sd_Control_pre = Control_pre_sd.Grand_RT_sd_pre;
sd_Control_mid = Control_mid_sd.Grand_RT_sd_mid;
sd_Control_post = Control_post_sd.Grand_RT_sd_post;

sd_Control = [sd_Control_pre', sd_Control_mid', sd_Control_post'];
sd_Retreat = [sd_Retreat_pre', sd_Retreat_mid', sd_Retreat_post'];

super_max_sd = max([max(sd_Control(:)), max(sd_Retreat(:))]);
super_min_sd = min([min(sd_Control(:)), min(sd_Retreat(:))]);

figure;
subplot(1,3,1);
errorbar(data_Control(1,:),sd_Control(1,:) );
hold on;
errorbar(data_Retreat(1,:),sd_Retreat(1,:));
legend({'Control','Retreat'},'Location','northeast');
title('Incongruent RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min-super_max_sd super_max+super_max_sd+10]);
xlim([0.8 3.2]);


subplot(1,3,2);
errorbar(data_Control(3,:), sd_Control(3,:));
hold on;
errorbar(data_Retreat(3,:), sd_Retreat(3,:));
legend({'Control','Retreat'},'Location','northeast');
title('Neutral RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min-super_max_sd super_max+super_max_sd+10]);
xlim([0.8 3.2]);


subplot(1,3,3);
errorbar(data_Control(5,:), sd_Control(5,:));
hold on;
errorbar(data_Retreat(5,:), sd_Retreat(5,:));
legend({'Control','Retreat'},'Location','northeast');
title('Congruent RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min-super_max_sd super_max+super_max_sd+10]);
xlim([0.8 3.2]);

%% do a bar plot: https://in.mathworks.com/matlabcentral/answers/379570-how-can-i-place-my-error-bar-in-separate-bar-center

data_incong = [data_Control_pre(1) data_Control_mid(1) data_Control_post(1); data_Retreat_pre(1) data_Retreat_mid(1) data_Retreat_post(1)]';
data_neu = [data_Control_pre(3) data_Control_mid(3) data_Control_post(3); data_Retreat_pre(3) data_Retreat_mid(3) data_Retreat_post(3)]';
data_cong = [data_Control_pre(5) data_Control_mid(5) data_Control_post(5); data_Retreat_pre(5) data_Retreat_mid(5) data_Retreat_post(5)]';

sd_incong = [sd_Control_pre(1) sd_Control_mid(1) sd_Control_post(1); sd_Retreat_pre(1) sd_Retreat_mid(1) sd_Retreat_post(1)]';
sd_neu = [sd_Control_pre(3) sd_Control_mid(3) sd_Control_post(3); sd_Retreat_pre(3) sd_Retreat_mid(3) sd_Retreat_post(3)]';
sd_cong = [sd_Control_pre(5) sd_Control_mid(5) sd_Control_post(5); sd_Retreat_pre(5) sd_Retreat_mid(5) sd_Retreat_post(5)]';

figure;
subplot(1,3,1);
hB=bar(data_incong);
X=[];
for i=1:length(hB)
  X=[X;hB(i).XData+hB(i).XOffset];
end
X = X.';     % rearrange X by column and show what this is...
hold on                                % so can add to bar plot
hEB=errorbar(X,data_incong, sd_incong,'.');  % add error bars
legend({'Control','Retreat'},'Location','northeast');
title('Incongruent RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min-super_max_sd super_max+super_max_sd+10]);
%xlim([0.8 3.2]);


subplot(1,3,2);
hB=bar(data_neu);
X=[];
for i=1:length(hB)
  X=[X;hB(i).XData+hB(i).XOffset];
end
X = X.';     % rearrange X by column and show what this is...
hold on                                % so can add to bar plot
hEB=errorbar(X,data_neu,sd_neu,'.');  % add error bars
legend({'Control','Retreat'},'Location','northeast');
title('Neutral RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min-super_max_sd super_max+super_max_sd+10]);
%xlim([0.8 3.2]);


subplot(1,3,3);
hB=bar(data_cong);
X=[];
for i=1:length(hB)
  X=[X;hB(i).XData+hB(i).XOffset];
end
X = X.';     % rearrange X by column and show what this is...
hold on                                % so can add to bar plot
hEB=errorbar(X,data_cong,sd_cong,'.');  % add error bars
legend({'Control','Retreat'},'Location','northeast');
title('Cong RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min-super_max_sd super_max+super_max_sd+10]);
%xlim([0.8 3.2]);
