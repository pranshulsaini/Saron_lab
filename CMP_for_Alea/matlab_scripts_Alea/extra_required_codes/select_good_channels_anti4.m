function [new_elp_file] = select_good_channels_anti4(full_elp_file, drop_list)

if ~isempty(drop_list) && strcmp(drop_list{1}, ' ')
    drop_list = [];
end


%Remember: drop_list is cell array of channels that are bad. 
fid = fopen(full_elp_file, 'r');

[pathstr, name, ext] = fileparts(full_elp_file);
new_elp_file = fullfile(pathstr, strcat('goodCh_', name,ext));
besa_elp_file = fullfile(pathstr, strcat('forBesa_', name,ext));
besa_ela_file = fullfile(pathstr, strcat('forBesa_', name,'.ela'));

fid2 = fopen(new_elp_file, 'w');
fid3 = fopen(besa_elp_file, 'w');
fid4 = fopen(besa_ela_file, 'w');

counter = 1; 
drop_list_num = [];
while(1)
    ln = fgets(fid);
    if ln == -1 
        break;
    end
    ln = strtrim(ln);
    if length(ln) < 3
        continue; 
    elseif length(ln) == 1   % this is to avoid blank lines
        continue;
    elseif strcmpi(strtrim(ln(1:3)), 'fid') == 1
        fprintf(fid3, '%s\n', ln);
        continue;
    elseif strcmpi(strtrim(ln(1:3)), 'REF') == 1
        fprintf(fid2, '%s\n', ln); 
        fprintf(fid3, '%s\n', ln); 
        fprintf(fid4, 'EEG REF_AVR\n');         
        continue;        
    end
    
    % Now checking for drop_list and writing new sfp
    match = 0;
    for i = 1:1:length(drop_list)
        num_chars=length(drop_list(i));
        if strcmp(strtrim(ln(1:num_chars)), drop_list(i)) == 1
            match = 1;
            drop_list_num = [drop_list_num counter]; %#ok<AGROW>
            break;
        end
    end
    
    % good channels, i.e. not found in the drop_list.
    if match == 0
        fprintf(fid2, '%s\n', ln);
        fprintf(fid3, '%s\n', ln);        
        if strcmpi(strtrim(ln(3)), '_') == 1
            fprintf(fid4, '%s\n', cell2mat(strcat('EEG',{' '},ln(1:2),{'_AVR'})));
        elseif strcmpi(strtrim(ln(4)), '_') == 1
            fprintf(fid4, '%s\n', cell2mat(strcat('EEG',{' '},ln(1:3),{'_AVR'})));
        else
            fprintf(fid4, '%s\n', cell2mat(strcat('EEG',{' '},ln(1:4),{'_AVR'})));
        end
    end
    counter = counter + 1;
end
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
