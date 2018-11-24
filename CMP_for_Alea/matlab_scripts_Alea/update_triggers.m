clc;
clear all;
%online link: https://www.mathworks.com/matlabcentral/answers/306876-how-do-i-read-only-a-specific-line-while-reading-a-text-file-in-matlab
%online link: https://www.mathworks.com/matlabcentral/answers/339342-how-to-read-a-letter-from-a-string-in-a-cell

sub_IDs = importdata('C:\Users\plsaini\Box Sync\Stroop\Subject_IDs.txt');
num_sub = size(sub_IDs,1);
cond = zeros(num_sub,1);

% 8 condition sequences. There are eight different orders of sequences and each has been assigned four columns in the excel file. The first column is
% the condition{-1,0,-1} for trials, the second one is the block number [1-16], the third one as block type {-1,0,1}, and the fourth one as blank col
conditions = importdata('C:\Users\plsaini\Box Sync\Stroop\conditions_sorted.xlsx');
conditions = cell2mat(struct2cell(conditions));

for i = 1:num_sub % 2nd file has incomplete event file. So a mismatch arises b/w log and event file
    i
    sub_IDs{i} = strrep(sub_IDs{i},'.bdf','');
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\Logfiles\',sub_IDs{i},'.txt'); % log file
    fid=fopen(filename);
    linenum = 7;      % This line number has condition written
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    cond_num = C{1}{1}(29);  % finding condition number
    cond(i) = str2num(cond_num);
    
    
    
    %link: https://www.mathworks.com/help/matlab/ref/textscan.html
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\pre_SOBI\',sub_IDs{i},'_aux_NoBad_AvgRef_OrigTrigs.evt');
    fid = fopen(filename,'r');
    MyText = textscan(fid,'%d %d %d %s %d','headerlines',1);   % will scan the whole text of the event file. First line will be assumed as heading
    fclose(fid);
    
    %log_file
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\Logfiles\',sub_IDs{i},'.txt'); 
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
    
    n_rows_stim = 0;
    ind = [];
    for p = last_31+1: n_rows_evt
        if ((MyText{1,3}(p)== 12) || (MyText{1,3}(p)== 5) || (MyText{1,3}(p)== 6) || (MyText{1,3}(p)== 20))  % only copy stimuli from event files
            
            % to check if the stimuli presented are same across logfile and event file
            if color2code_stim(log_file{1,3}(s)) ~= MyText{1,3}(p)
                disp('ERROR: Stimuli are not matched between log and event file. Check!!!!!');
                ind(end+1) = s; % I have to remove this condition from the stimuli sequences later on
                s = s +1;
                
                if color2code_stim(log_file{1,3}(s)) == MyText{1,3}(p)
                    disp('The problem is fixed')
                else
                    disp('The problem not solved');
                end
   
            end
            
    
            col1(q) = MyText{1,1}(p);  % stimulus time stamp
            col1(r) = col1(q) + 1000* log_file{1,5}(s);  % bcz log file reaction time is in ms
            
            col2(q) = MyText{1,2}(p);  % '1'
            col2(r) = col2(q); % '1'
            
            col3(q) = MyText{1,3}(p);    % stimulus
            col3(r) = color2code(log_file{1,4}(s));  % response
            if(col3(r) == 500) % miss trial
                outcome(s) = 500;  % otherwise outcome was just correct/incorrect
                col1(r) = col1(r) + 1000* max(log_file{1,5}(:)); % because the reaction time for miss trial was coded as 0
            end
            
            col4(q) = MyText{1,4}{p}; % 'Trig'
            col4(r) = col4(q);  % 'Trig'
            
            col5(q) = MyText{1,5}(p); % stimulus
            col5(r) = col3(r); % response
            
             
            
            q = q+2;  % will write only in odd rows
            r = r+2; % will write only in even rows 
            s = s+1; 
            n_rows_stim = n_rows_stim +1; % to find out the number of stimuli found in the original event file (usually they would be 256 after the loop ends)
           
        end
    end
   
    if s ~= 257
        disp('The number of stimuli in event file are not 256. This is erroneous');
    end
    
    col4 = cell2mat(col4);
    col5 = num2str(col5);
    
    % storing old triggers
    old_triggers = col3;
    
    
    % writing the corrected original trigger file where we have extracted
    % the missing responses from the log file.
    name = strcat('C:\Users\plsaini\Box Sync\Stroop\pre_SOBI\', sub_IDs{i},'_orig_triggers');
    save(name,'old_triggers');
    
    n_rows = n_rows_stim*2; % for every trial, we have two triggers
    
    % updating triggers. 
    cond_seq = (conditions(:,4*cond(i)-3) + conditions(:,4*cond(i)-1) + 3)*1000; % the sum of 1st two terms lead to -2(incongruent), -1 (neutral in incongruent), 0 (neutral in neutral)
    %, 1 (neutral in congruent), 2 congruent. After adding 3, now it will be 1000,2000,3000,4000,5000 
    
    %Now we will encode trial numbers in this sequence
    cond_seq = cond_seq + repmat(1:128,1,2)';
    
    %deleting those conditions whose corresponding stimuli were not in event file. This happens rarely but some files are screwed up.
    cond_seq(ind) = [];
    
    %deleting responses corresponding to non-existing stimuli in the original event file
    outcome(ind) = [];
    
    col3(1:2:n_rows) = cond_seq;  % updating values for stimulus rows
    col3(2:2:n_rows) = outcome;  % it accounts for none values as well
    col5 = string(col5); % necessary because char arrays are two dimensional which causes problems in copying values from outcome and cond_seq
    col5(1:2:n_rows) = string(cond_seq); % updating values for stimulus rows
    col5(2:2:n_rows) = string(outcome);
    col5 = char(col5); % reconversion otherwise there was a problem in print file
    
    

    
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\pre_SOBI\',sub_IDs{i},'_aux_NoBad_AvgRef_UpdatedTrigs.evt');
    fid = fopen(filename,'w');
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
 
    for j = 1:n_rows
        fprintf(fid, '%d\t%d\t%d\t%s\n', col1(j), col2(j), col3(j), [col4(j,:), col5(j,:)]);
    end
    
    fclose(fid);
    

    
end




