%% This file will average reaction time data across subjects for each group, condition, and assessment

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
RT_avg_pre = zeros(round(num_sub/3),2); % 5 conditions
RT_avg_mid = zeros(round(num_sub/3),2); % 5 conditions
RT_avg_post = zeros(round(num_sub/3),2); % 5 conditions

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
    
    order = rem(cond(i),2);  % this will help in finding out the trial numbers corresponding to the neutral blocks
    
    if order == 0
        neu_trial = [1:4:16 34:4:48 67:4:80 100:4:112 129:4:144 162:4:176 195:4:208 228:4:240 ];
    elseif order == 1
        neu_trial = [17:4:32 50:4:64 83:4:96 116:4:128 145:4:160 178:4:192 211:4:224 244:4:256];
    else
        disp('unvalid sequence order')
    end
    
    
    % RT averageing is being done for pre and only for correct trials
    if (sub_IDs{i}(8) == '1')
        count_pre = count_pre + 1;
        RT_neu = [];
        RT_neu32 = [];

        for j = 1:256
            if ((cond_seq(j) == 3000) && (outcome(j) == 400))
                RT_neu(end + 1) = log_file{1,5}(j);
            end            
        end
        
        for j = 1:32
            if (outcome(neu_trial(j)) == 400)
                RT_neu32(end + 1) = log_file{1,5}(neu_trial(j));
            end
        end
                
        
    RT_avg_pre(count_pre,:) =[median(RT_neu), median(RT_neu32)];
    counter = 1;
    end
    
    
    % RT averageing is being done for mid and only for correct trials
    if (sub_IDs{i}(8) == '2')
        count_mid = count_mid + 1;
        RT_neu = [];
        RT_neu32 = [];

        for j = 1:256
            if ((cond_seq(j) == 3000) && (outcome(j) == 400))
                RT_neu(end + 1) = log_file{1,5}(j);
            end            
        end
        
        for j = 1:32
            if (outcome(neu_trial(j)) == 400)
                RT_neu32(end + 1) = log_file{1,5}(neu_trial(j));
            end
        end
                
        
    RT_avg_mid(count_mid,:) =[median(RT_neu), median(RT_neu32)];
    counter = 1;

    end
    
    
    % RT averageing is being done for post and only for correct trials
    if (sub_IDs{i}(8) == '3')
        count_post = count_post + 1;
        RT_neu = [];
        RT_neu32 = [];

        for j = 1:256
            if ((cond_seq(j) == 3000) && (outcome(j) == 400))
                RT_neu(end + 1) = log_file{1,5}(j);
            end            
        end
        
        for j = 1:32
            if (outcome(neu_trial(j)) == 400)
                RT_neu32(end + 1) = log_file{1,5}(neu_trial(j));
            end
        end
                
        
    RT_avg_post(count_post,:) =[median(RT_neu), median(RT_neu32)];
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

save(strcat(group,'_RT_avg_pre_median_neu_32_128'),'RT_avg_pre');
save(strcat(group,'_RT_avg_mid_median_neu_32_128'),'RT_avg_mid');
save(strcat(group,'_RT_avg_post_median_neu_32_128'),'RT_avg_post');

% averaging across subjects as well
Grand_RT_avg_pre = median(RT_avg_pre);
Grand_RT_avg_mid = median(RT_avg_mid);
Grand_RT_avg_post = median(RT_avg_post);

save(strcat(group,'_Grand_RT_avg_pre_median_neu_32_128'),'Grand_RT_avg_pre');
save(strcat(group,'_Grand_RT_avg_mid_median_neu_32_128'),'Grand_RT_avg_mid');
save(strcat(group,'_Grand_RT_avg_post_median_neu_32_128'),'Grand_RT_avg_post');


%% assessment wise plotting
figure;
subplot(3,1,1);
plot(RT_avg_pre);
%legend({'incong','incong-neu','neu','cong-neu','cong'},'Location','southeast')
legend({'neu128','neu32'},'Location','southeast')
title(strcat(group,' RT at pre'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

subplot(3,1,2);
plot(RT_avg_mid);
%legend({'incong','incong-neu','neu','cong-neu','cong'},'Location','southeast')
legend({'neu128','neu32'},'Location','southeast')
title(strcat(group,' RT at mid'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

subplot(3,1,3)
plot(RT_avg_post);
%legend({'incong','incong-neu','neu','cong-neu','cong'},'Location','southeast')
legend({'neu128','neu32'},'Location','southeast')
title(strcat(group,' RT at post'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

%% condition-wise plotting
figure;
subplot(2,1,1);
plot(RT_avg_pre(:,1));
hold on;
plot(RT_avg_mid(:,1));
plot(RT_avg_post(:,1));
legend({'pre','mid','post'},'Location','southeast')
title(strcat(group,' RT: neutral128'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);

subplot(2,1,2);
plot(RT_avg_pre(:,2));
hold on;
plot(RT_avg_mid(:,2));
plot(RT_avg_post(:,2));
legend({'pre','mid','post'},'Location','southeast')
title(strcat(group,' RT: neutral32'));
xlabel('sub-number');
ylabel('RT (ms)');
ylim([super_min super_max]);


%% plotting grand RT averages
Retreat_pre = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_avg_pre_median_neu_32_128.mat');
Retreat_mid = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_avg_mid_median_neu_32_128.mat');
Retreat_post = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\retreat\Retreat_Grand_RT_avg_post_median_neu_32_128.mat');
Control_pre = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_avg_pre_median_neu_32_128.mat');
Control_mid = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_avg_mid_median_neu_32_128.mat');
Control_post = load('C:\Users\plsaini\Box Sync\Stroop\Behaviour\RT_results\control\Control_Grand_RT_avg_post_median_neu_32_128.mat');

data_Retreat_pre = Retreat_pre.Grand_RT_avg_pre;
data_Retreat_mid = Retreat_mid.Grand_RT_avg_mid;
data_Retreat_post = Retreat_post.Grand_RT_avg_post;
data_Control_pre = Control_pre.Grand_RT_avg_pre;
data_Control_mid = Control_mid.Grand_RT_avg_mid;
data_Control_post = Control_post.Grand_RT_avg_post;

data_Control = [data_Control_pre', data_Control_mid', data_Control_post'];
data_Retreat = [data_Retreat_pre', data_Retreat_mid', data_Retreat_post'];

super_max = max([max(data_Control(:)), max(data_Retreat(:))]);
super_min = min([min(data_Control(:)), min(data_Retreat(:))]);

figure;
subplot(1,2,1);
plot(data_Control(1,:));
hold on;
plot(data_Control(2,:));
legend({'Neu128','Neu32'},'Location','northeast');
title('Control RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min super_max]);

subplot(1,2,2);
plot(data_Retreat(1,:));
hold on;
plot(data_Retreat(2,:));
legend({'Neu128','Neu32'},'Location','northeast');
title('Retreat RT Average');
xlabel('Assessment(Pre, Mid, Post)')
ylabel('RT (ms)');
ylim([super_min super_max]);


