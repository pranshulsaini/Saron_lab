function MSI_ANTI_Permutation1(xlsfile,percent_data,freq_bin_calc,vol_cond) 
% updtaed 02/01/2012
%function MSI_ANTI_Permutation1(xlsfile,percent_data,freq_bin_calc,vol_cond) 
  
  [num,txt1,txt2]=xlsread(xlsfile); % txt1 contains string and txt2 contains numeric values
  [subjects_no, para_no]=size(txt1);
  subjects_no=subjects_no-1; % considering the header of the xlsfile
  counter_set=0;
   
  for i=2:subjects_no+1
    flag=cell2num(txt2(i,19));
    if flag==1             % flag means the number of trials in this dataset is greater than minimum acceptable threshold
        counter_set=counter_set+1;
        dataset=cell2mat(strcat(txt1(i,3),'\',txt1(i,8)));
        traillength=cell2num(txt2(i,12));
        evtfile=cell2num(strcat(txt1(i,3),'\',txt1(i,9)));
        sfpfile=cell2mat(strcat(txt1(i,3),'\',txt1(i,10)));
        Fs=cell2num(txt2(i,13));
        prestim=cell2num(txt2(i,14));
        stim_trig=cell2num(txt2(i,4));
        res_trig=cell2num(txt2(i,5));
        elec_num=cell2num(txt2(i,15));
        working_folder=cell2mat(txt1(i,3));
        EEG_datasetname=strcat(cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))))
        [trialdata1 cfg timelock]= MSI_ANTI_permutation_all_steps1(dataset, traillength,evtfile, sfpfile,Fs,prestim,stim_trig,res_trig,elec_num,working_folder,EEG_datasetname);
                 
  end 
  end

     