% This file will extract approx 32 trials out of 128 pure neutral trials we have

% %Steps
% 1. Read the sequence order number from the log file of the subject
% 2. Even order numbers start with neutral blocks where as odd order numbers start with congruent/incongruent blocks. This information tells us the trial numbers of neutral blocks
% 3. Go the beginning of each neutal block extract the information of all the trials I want.
% 4. Copy the above information into a new file. The averages will be created only for these triggers now
%%

clc;
%online link: https://www.mathworks.com/matlabcentral/answers/306876-how-do-i-read-only-a-specific-line-while-reading-a-text-file-in-matlab
%online link: https://www.mathworks.com/matlabcentral/answers/339342-how-to-read-a-letter-from-a-string-in-a-cell

sub_IDs = importdata('C:\Users\plsaini\Box Sync\Stroop\Subject_IDs.txt');
num_sub = size(sub_IDs,1);
cond = zeros(num_sub,1);

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
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\STR_SOBI_Output\',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI\',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_Sources_STR\',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_recon\','recon_',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_withRespTrigs.evt');
    fid = fopen(filename,'r');
    MyText = textscan(fid,'%d %d %d %s','headerlines',1);   % will scan the whole text of the event file. First line will be assumed as heading
    fclose(fid);
    
    
    order = rem(cond(i),2);  % this will help in finding out the trial numbers corresponding to the neutral blocks
    
    if order == 0
        neu_trial = [1:4:16 34:4:48 67:4:80 100:4:112];
    elseif order == 1
        neu_trial = [17:4:32 50:4:64 83:4:96 116:4:128];
    else
        disp('unvalid sequence order')
    end
    
    for k = 1:128  % loop is necessary because 
        index = find(rem(MyText{1,3}(:),1000) == k); % It will tell when the next set of [1-128] trials begin
        if length(index) == 1 && index>100  % 100 is chosen as an approximate number which would make sure that we are in the second set
            first_set_end = index - 1;
            second_set_start = index;
            second_set_end = length(MyText{1,3}(:));
            break;
        end
            
        if length(index) == 2
            first_set_end = index(2) - 1;
            second_set_start = index(2);
            second_set_end = length(MyText{1,3}(:));
            break;
        end
    end
    

    start = 1;
    
    filename = strcat('C:\Users\plsaini\Box Sync\Stroop\STR_SOBI_Output\',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI\',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_Sources_STR\',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_recon\','recon_',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_neutralwithRespTrigs.evt');
    fid = fopen(filename,'w');
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    
    prvs_tmp = 0; % value of pervious timestamp. WIll be used to keep check of some errors
    
    for i = 1:16
        for j = start: first_set_end
            if (rem(MyText{1,3}(j),1000) == neu_trial(i))
                fprintf(fid, '%d\t%d\t%d\t%s\n', MyText{1,1}(j), MyText{1,2}(j), MyText{1,3}(j), cell2mat(MyText{1,4}(j))); % for stimulus trigger
                fprintf(fid, '%d\t%d\t%d\t%s\n', MyText{1,1}(j+1), MyText{1,2}(j+1), MyText{1,3}(j+1), cell2mat(MyText{1,4}(j+1))); % for response trigger
                diff = MyText{1,1}(j) - prvs_tmp;
                if diff < 0
                    disp('There is an error in the order of copying trials')
                end
                prvs_tmp = MyText{1,1}(j+1);
                start = j+2;
                break;
            end
        end
    end
    
    for i = 1:16
        for j = start: second_set_end
            if (rem(MyText{1,3}(j),1000) == neu_trial(i))
                fprintf(fid, '%d\t%d\t%d\t%s\n', MyText{1,1}(j), MyText{1,2}(j), MyText{1,3}(j), cell2mat(MyText{1,4}(j))); % for stimulus trigger
                fprintf(fid, '%d\t%d\t%d\t%s\n', MyText{1,1}(j+1), MyText{1,2}(j+1), MyText{1,3}(j+1), cell2mat(MyText{1,4}(j+1))); % for response trigger
                diff = MyText{1,1}(j) - prvs_tmp;
                if diff < 0
                    disp('There is an error in the order of copying trials')
                end
                prvs_tmp = MyText{1,1}(j+1);
                start = j+2;
                break;
            end
        end
    end
    
    
    fclose(fid);

end




