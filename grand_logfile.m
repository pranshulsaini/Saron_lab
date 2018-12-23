% Creating an excel file as asked by Brandon

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
    sub_IDs{i} = strrep(sub_IDs{i},'.txt','');
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\Logfiles\',sub_IDs{i},'.txt'); % log file
    fid=fopen(filename);
    linenum = 7;      % This line number has condition written
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    cond_num = C{1}{1}(29);  % finding condition number
    cond(i) = str2num(cond_num);
    fclose(fid);

    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\Logfiles\',sub_IDs{i},'.txt'); 
    fid = fopen(filename);
    log_file = textscan(fid,'%d %s %s %s %d %d','headerlines',98);
    fclose(fid);
    
    col_header = {'Filename','Subject','Group','Retreat','Testing','Trial','EEG_File','Block_Num','Block_Type','Trial_Type','Word','Color',	'Response', 'Correct','Miss','RT','ITI'};
    
    filenames = strings(256,1);
    filenames(:) = sub_IDs{i};
    
    subject = zeros(256,1);
    subject(:) = str2num(sub_IDs{i}(4:6));
    
    group = strings(256,1);
    group(:) =   'retreat';
    
    retreat = ones(256,1);  % retreat 1

    testing = strings(256,1);
    if sub_IDs{i}(8) == '1' 
        val = 'pre';
    elseif sub_IDs{i}(8) == '2'
        val = 'mid';
    elseif sub_IDs{i}(8) == '3'
        val = 'post';
    else
        display('invalid value for assessment');
    end
    testing(:) = val;
    
    trial = [1:256]';
    
    eeg_file = ones(256,1); % In this code, I am using only the logfiles which correspond to accepted EEG files
    
    block_num = fix([0:255]'/16) + 1;   % fix gives out only the integer part
    
    block_type = conditions(:,4*cond(i)-1);
    
    trial_type = conditions(:,4*cond(i)-3);
    
    word = strings(256,1);
    word(:) = log_file{1,2}(:);
    
    color =  strings(256,1);
    color(:) = log_file{1,3}(:);
    
    color_len = cellfun('length',color);
    
    response =  strings(256,1);
    response(:) = log_file{1,4}(:);  
 
    resp_len = cellfun('length', response);
    
    correct = double(color_len == resp_len);
    correct(response == 'none') = 0;
    
    miss = zeros(256,1);
    miss(:) = response == 'none';
    
    rt = zeros(256,1);
    rt(:) = log_file{1,5}(:); 
    
    iti = zeros(256,1);
    iti(:) = log_file{1,6}(:);
       
    
    filename = 'C:\Users\plsaini\Box Sync\Stroop\STR_logfile_grand.xlsx';
    start = (i-1)*256 + 2;
    endd = i*256 + 1; 
    xlswrite(filename,col_header,'Sheet1','A1');     %Write column header
    
    xlswrite(filename,filenames,'Sheet1',strcat('A',num2str(start),':','A',num2str(endd)));   % write data
    xlswrite(filename,subject,'Sheet1',strcat('B',num2str(start),':','B',num2str(endd)));   % write data
    xlswrite(filename,group,'Sheet1',strcat('C',num2str(start),':','C',num2str(endd)));   % write data
    xlswrite(filename,retreat,'Sheet1',strcat('D',num2str(start),':','D',num2str(endd)));   % write data
    xlswrite(filename,testing,'Sheet1',strcat('E',num2str(start),':','E',num2str(endd)));   % write data
    xlswrite(filename,trial,'Sheet1',strcat('F',num2str(start),':','F',num2str(endd)));   % write data
    xlswrite(filename,eeg_file,'Sheet1',strcat('G',num2str(start),':','G',num2str(endd)));   % write data
    xlswrite(filename,block_num,'Sheet1',strcat('H',num2str(start),':','H',num2str(endd)));   % write data
    xlswrite(filename,block_type,'Sheet1',strcat('I',num2str(start),':','I',num2str(endd)));   % write data
    xlswrite(filename,trial_type,'Sheet1',strcat('J',num2str(start),':','J',num2str(endd)));   % write data
    xlswrite(filename,word,'Sheet1',strcat('K',num2str(start),':','K',num2str(endd)));   % write data
    xlswrite(filename,color,'Sheet1',strcat('L',num2str(start),':','L',num2str(endd)));   % write data
    xlswrite(filename,response,'Sheet1',strcat('M',num2str(start),':','M',num2str(endd)));   % write data
    xlswrite(filename,correct,'Sheet1',strcat('N',num2str(start),':','N',num2str(endd)));   % write data
    xlswrite(filename,miss,'Sheet1',strcat('O',num2str(start),':','O',num2str(endd)));   % write data
    xlswrite(filename,rt,'Sheet1',strcat('P',num2str(start),':','P',num2str(endd)));   % write data
    xlswrite(filename,iti,'Sheet1',strcat('Q',num2str(start),':','Q',num2str(endd)));   % write data

end

