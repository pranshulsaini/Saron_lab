clc;
clear all;

mainpath = 'C:\Users\plsaini\Box Sync\Stroop\STR_SOBI_Output\';
sub_IDs = importdata('C:\Users\plsaini\Box Sync\Stroop\Subject_IDs.txt'); % It has the names of all the subjects (in a column) who need to be processed, for example: STR00111
num_sub = size(sub_IDs,1);

for i = 1:num_sub % 2nd file has incomplete event file. So a mismatch arises b/w log and event file
    i
    sub_IDs{i} = strrep(sub_IDs{i},'.bdf','');
    localpath = strcat(sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI\', sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_Sources_STR\',sub_IDs{i}, '_aux_NoBad_AvgRef_forSOBI_recon\');
    
    filename = strcat(mainpath,localpath,'recon_',sub_IDs{i},'_aux_NoBad_AvgRef_forSOBI_av-export.mul'); % the address and filename for .mul file

    data = importdata(filename);
    data = data.data;

    n_chan = size(data,2);
    n_iter = round(size(data,1)/4096);

    for j = 1:n_iter
        start = 1 + (j-1)*4096;
        endd =  j*4096;
        data_for_eph = data(start:endd,:);

        filename = strcat(mainpath,localpath,sub_IDs{i},'_cond_', num2str(j),'.eph');

        row1 = [n_chan, 4096, 2048];  % 1st element are the number of channels, 2nd is the total time points, and 3rd one is sampling rate

        dlmwrite(filename,row1,'delimiter','\t');
        dlmwrite(filename,data_for_eph','-append','delimiter','\t');
    end
end
