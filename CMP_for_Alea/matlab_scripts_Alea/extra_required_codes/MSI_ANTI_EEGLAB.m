function MSI_ANTI_EEGLAB(xlsfile,percent_data,freq_bin) 
%function [ALLEEG,EEG,CURRENTSET]= MSI_ANTI_EEGLAB(xlsfile,percent_data,freq_bin) 

  [num,txt1,txt2]=xlsread(xlsfile); % txt1 contains string and txt2 contains numeric values
  [subjects_no, para_no]=size(txt1);
  subjects_no=subjects_no-1; % considering the header of the xlsfile
   
  for i=2:subjects_no+1
    flag=cell2num(txt2(i,19));
    if flag==1             % flag means the number of trials in this dataset is greater than minimum acceptable threshold
        
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
        EEG_datasetname=strcat(cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
        [trialdata1]= MSI_ANTI_all_steps_fieldtrip6(dataset, traillength,evtfile, sfpfile,Fs,prestim,stim_trig,res_trig,elec_num,working_folder,EEG_datasetname);
         % writing triiger event file over the xls file
        range = sprintf('R%i', i);
        output={strcat(EEG_datasetname,'_Stim_Trigger.evt')};
        stim_trig_file=strcat(EEG_datasetname,'_Stim_Trigger.evt');
        xlswrite(xlsfile,output, 1, range);
       
    
        [ALLEEG EEG] = fieldtrip2eeglab6(trialdata1);
        EEG = pop_select(EEG, 'nochannel', cell2num(txt2(i,15))+1);
        EEG = pop_editset( EEG, 'setname', EEG_datasetname,'xmin',-prestim/Fs);
        EEG.chanlocs = readlocs( strcat(working_folder,'\' ,cell2mat(txt1(i,11))), 'filetype', 'sfp'); 
        EEG = pop_importepoch( EEG, strcat(working_folder,'\' ,stim_trig_file), {'Time' ,'Trig'}, 'latencyfields',{'Time'}, 'typefield', 'Trig', 'timeunit',1, 'headerlines',1);
        EEG = pop_importepoch( EEG, strcat(working_folder,'\' ,cell2mat(txt1(i,17))), {'Response_Time'}, 'latencyfields',{'Response_Time'},  'timeunit',1e-6, 'headerlines',1);
    
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,i);
     
        EEG = pop_saveset( EEG, 'filepath', working_folder,'filename',EEG_datasetname);
        [status]=coherence_raster6(EEG,EEG_datasetname, working_folder, percent_data,freq_bin);
        if status==1
           display(strcat('Coherence Maps have been generated for dataset....',EEG_datasetname));
        end 
        
    end
    
  end
end
