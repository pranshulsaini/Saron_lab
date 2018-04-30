clc;
clear all;
%online link: https://www.mathworks.com/matlabcentral/answers/306876-how-do-i-read-only-a-specific-line-while-reading-a-text-file-in-matlab
%online link: https://www.mathworks.com/matlabcentral/answers/339342-how-to-read-a-letter-from-a-string-in-a-cell

sub_IDs = importdata('D:\Stroop\Temp\Subject_IDs.txt');
num_sub = size(sub_IDs,1);
cond = zeros(num_sub,1);

% 8 condition sequences
conditions = importdata('D:\Stroop\Temp\conditions_sorted.xlsx');
conditions = cell2mat(struct2cell(conditions));

for i = 1:num_sub % 2nd file has incomplete event file. So a mismatch arises b/w log and event file
    i
    filename = strcat('D:\Stroop\Temp\log_files\',sub_IDs{i},'.txt'); % log file
    fid=fopen(filename);
    linenum = 7;      % This line number has condition written
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    cond_num = C{1}{1}(29);  % finding condition number
    cond(i) = str2num(cond_num);
    
    % finding 
    
    
    %link: https://www.mathworks.com/help/matlab/ref/textscan.html
    filename = strcat('D:\Stroop\Temp\pre_sobi\',sub_IDs{i},'_aux_NoBad_AvgRef_OrigTrigs.evt');
    fid = fopen(filename,'r');
    MyText = textscan(fid,'%d %d %d %s %d','headerlines',1);   % will scan the whole text of the event file. First line will be assumed as heading
    fclose(fid);
    
    %log_file
    filename = strcat('D:\Stroop\Temp\log_files\',sub_IDs{i},'.txt'); 
    fid = fopen(filename);
    log_file = textscan(fid,'%d %s %s %s %d %d','headerlines',98);
    stim = log_file{1,3}(:);
    stim_len = cellfun('length',stim);
    resp = log_file{1,4}(:);
    resp_len = cellfun('length', resp);
    outcome = stim_len == resp_len;
    outcome = outcome + 399; % 400 for correct, 399 for incorrect
       
    
    
    % finding the starting point for copying
    index = find(MyText{1,3}(:) == 31);  % finding index of all 31 triggers 
    last_31 = index(end);  % index of last 31
    
    
    %defining the columns where data has to be copied from event files
    n_rows_log = length(stim); 
    col1 = zeros(2*n_rows_log,1); % two times because the event files has stimuli and responses separately listed
    col2 = zeros(2*n_rows_log,1);
    col3 = zeros(2*n_rows_log,1);
    %col4 = char(zeros(n_rows_log,5));  
    col4 = strings(2*n_rows_log,1);
    col5 = zeros(2*n_rows_log,1);
    
    %copy data to the columns 
    n_rows_evt = length(MyText{1,3}(:));
    q= 1; % to write stimuli in column files
    r = 2; % to write responses in column files
    s= 1;  % to extract responses from logfile
    for p = last_31+1: n_rows_evt
        if ((MyText{1,3}(p)== 12) || (MyText{1,3}(p)== 5) || (MyText{1,3}(p)== 6) || (MyText{1,3}(p)== 20))  % only copy stimuli from event files
            col1(q) = MyText{1,1}(p);  % stimulus time stamp
            col1(r) = col1(q) + 1000* log_file{1,5}(s);  % bcz log file reaction time is in ms
            
            col2(q) = MyText{1,2}(p);  % '1'
            col2(r) = col2(q); % '1'
            
            col3(q) = MyText{1,3}(p);    % stimulus
            col3(r) = color2code(log_file{1,4}(s));  % response
            if(col3(r) == 500)
                outcome(s) = 500;  % otherwise outcome was just correct/incorrect
            end
            
            col4(q) = MyText{1,4}{p}; % 'Trig'
            col4(r) = col4(q);  % 'Trig'
            
            col5(q) = MyText{1,5}(p); % stimulus
            col5(r) = col3(r); % response
            
             
            
            q = q+2;  % will write only in odd rows
            r = r+2; % will write only in even rows 
            s = s+1; 
        end
    end
   
    col4 = cell2mat(col4);
    col5 = num2str(col5);
    
    % storing old triggers
    old_triggers = col3;
    
    % updating triggers
    cond_seq = conditions(:,4*cond(i)-3) + conditions(:,4*cond(i)-1) + 203; % the sum of 1st two terms lead to -2(incongruent), -1 (neutral in incongruent), 0 (neutral in neutral), 1 (neutral in congruent), 2 congruent 
    
    col3(1:2:end) = cond_seq;  % updating values for stimulus rows
    col3(2:2:end) = outcome;  % it accounts for none values as well
    col5 = string(col5); % necessary because char arrays are two dimensional which causes problems in copying values from outcome and cond_seq
    col5(1:2:end) = string(cond_seq); % updating values for stimulus rows
    col5(2:2:end) = string(outcome);
    col5 = char(col5); % reconversion otherwise there was a problem in print file
    
    
    % adding trial triggers
    n_rows = n_rows_log*4;  % for every trial, we have 4 triggers: block number, trial number, condition, response
    
    col1_trial = int32(zeros(n_rows, 1));
    col2_trial = int32(zeros(n_rows, 1));
    col3_trial = int32(zeros(n_rows, 1));
    col4_trial = char(zeros(n_rows, 5));
    col5_trial = char(zeros(n_rows, 3));

    %copying values to new column vectors
    for n = 1:256
        j= 4*n -1;
        k = 2*n-1;
        col1_trial(j:j+1) = col1(k:k+1);
        col2_trial(j:j+1) = col2(k:k+1);
        col3_trial(j:j+1) = col3(k:k+1);
        col4_trial(j:j+1,:) = col4(k:k+1,:);
        col5_trial(j:j+1,:) = col5(k:k+1,:);
    end
    
    
    %adding block number triggers
    col1_trial(1:4:end) = col1(1:2:end);  % stimulus timestamp values
    col2_trial(:) = 1;
    col3_trial(1:4:end,:) = conditions(:,4*cond(i)-2);  % block number
    col4_trial =  string(col4_trial);
    col4_trial(:) = 'Trig.';
    col4_trial =  char(col4_trial);
    col5_trial = string(col5_trial);
    col5_trial(1:4:end) = string(conditions(:,4*cond(i)-2));  % block number
    col5_trial = char(col5_trial);
    
    
    
    %adding trial number triggers
    col1_trial(2:4:end) = col1(1:2:end);  % stimulus timestamp values
    col2_trial(:) = 1; 
    col3_trial(2:4:end,:) = repmat(101:116,1,16)' ; 
    col4_trial =  string(col4_trial);
    col4_trial(:) = 'Trig.';
    col4_trial =  char(col4_trial);
    col5_trial = string(col5_trial);
    col5_trial(2:4:end) = string(repmat(1:16,1,16)');
    col5_trial = char(col5_trial);
    
    
    filename = strcat(sub_IDs{i},'_aux_NoBad_AvgRef_UpdatedTrigs.evt');
    fid = fopen(filename,'w');
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');

    for j = 1:n_rows
        fprintf(fid, '%d\t%d\t%d\t%s\n', col1_trial(j), col2_trial(j), col3_trial(j), [col4_trial(j,:), col5_trial(j,:)]);
    end
    
    fclose(fid);
    
    name = strcat(sub_IDs{i},'_old_triggers');
    save(name,'old_triggers');
    
end




