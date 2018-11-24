function MSI_ANTI_EEGLAB2(xlsfile,percent_data,freq_bin_calc) 
%function [ALLEEG,EEG,CURRENTSET]= MSI_ANTI_EEGLAB(xlsfile,percent_data,freq_bin) 
  ALLTRAILS = eeg_emptyset();
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
       
        if counter_set==1
            [ALLTRAILS]=EEG;
        else 
            [ALLTRAILS]=pop_mergeset(ALLTRAILS,EEG);

        end 
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,i);
        EEG = pop_saveset( EEG, 'filepath', working_folder,'filename',EEG_datasetname);
     end
  end

if freq_bin_calc==1
   [alpha_range]=IAF_Calc6(ALLTRAILS);
   fid = fopen(strcat(working_folder,'\',EEG_datasetname,'_Freq_IAF.log'),'w');
   if fid == -1
      fprintf(1,'Error creating New Event File  \n');
   end
   fprintf(fid,'Alpha_F1\tAlpha_F2\tAlpha_F3\tAlpha_F4\n');
   fprintf(fid,'%d\t%d\t%d\t%d\n', alpha_range(1) ,alpha_range(2) ,alpha_range(3), alpha_range(4));
   fclose(fid);
end
freq_bin=[alpha_range(1) alpha_range(4)];

  
  for i=2:subjects_no+1
    flag=cell2num(txt2(i,19));
    if flag==1             % flag means the number of trials in this dataset is greater than minimum acceptable threshold
       working_folder=cell2mat(txt1(i,3));
       EEG_datasetname=strcat(cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
       [EEG]=pop_loadset(strcat(EEG_datasetname,'.set'), working_folder);
       [status]=coherence_raster6(EEG,EEG_datasetname, working_folder, percent_data,freq_bin);
        if status==1
           display(strcat('Coherence Maps have been generated for dataset....',EEG_datasetname));
        end 
    end
  end
  
end
