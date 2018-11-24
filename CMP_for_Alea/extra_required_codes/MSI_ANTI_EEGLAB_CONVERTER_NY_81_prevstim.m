function MSI_ANTI_EEGLAB_CONVERTER_NY_81_prevstim(xlsfile,with_res_trig) 
% updtaed 12/21/2012
%function [ALLEEG,EEG,CURRENTSET]= MSI_ANTI_EEGLAB(xlsfile,percent_data,freq_bin) 
  ALLTRAILS = eeg_emptyset(); % for concatenating all connditions for IAF estimation
  [num,txt1,txt2]=xlsread(xlsfile); % txt1 contains string and txt2 contains numeric values
  [subjects_no, para_no]=size(txt1);
  subjects_no=subjects_no-1; % considering the header of the xlsfile
  counter_set=0;
  close all; 
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
        prev_stim=cell2num(txt2(i,23));
        working_folder=cell2mat(txt1(i,3));
        EEG_datasetname=strcat(cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))),'_PrevStim_',mat2str(cell2mat(txt2(i,23))));
        if prev_stim>0
           [trialdata1]= MSI_ANTI_all_steps_fieldtrip6_NY_prevstim(dataset, traillength,evtfile, sfpfile,Fs,prestim,stim_trig,res_trig,elec_num,working_folder,EEG_datasetname,with_res_trig,prev_stim);
        else
           [trialdata1]= MSI_ANTI_all_steps_fieldtrip6_NY(dataset, traillength,evtfile, sfpfile,Fs,prestim,stim_trig,res_trig,elec_num,working_folder,EEG_datasetname,with_res_trig);
        end 
         % writing triiger event file over the xls file
        range = sprintf('R%i', i);
        output={strcat(EEG_datasetname,'_Stim_Trigger.evt')};
        stim_trig_file=strcat(EEG_datasetname,'_Stim_Trigger.evt');
        xlswrite(xlsfile,output, 1, range);
       
        [ALLEEG EEG] = fieldtrip2eeglab6(trialdata1);
        EEG = pop_select(EEG, 'nochannel', cell2num(txt2(i,15))+1); % removing EDF info Chn
        EEG = pop_editset( EEG, 'setname', EEG_datasetname,'xmin',-prestim/Fs);
        EEG.chanlocs = readlocs( strcat(working_folder,'\' ,cell2mat(txt1(i,11))), 'filetype', 'sfp'); 
        EEG = pop_importepoch( EEG, strcat(working_folder,'\' ,stim_trig_file), {'Time' ,'Trig'}, 'latencyfields',{'Time'}, 'typefield', 'Trig', 'timeunit',1, 'headerlines',1);
        if with_res_trig==1
            EEG = pop_importepoch( EEG, strcat(working_folder,'\' ,cell2mat(txt1(i,17))), {'Response_Time'}, 'latencyfields',{'Response_Time'},  'timeunit',1e-6, 'headerlines',1);
        end 
        
        if counter_set==1
            [ALLTRAILS]=EEG;
        else 
            [ALLTRAILS]=pop_mergeset(ALLTRAILS,EEG);

        end
         EEG.chanlocs = readlocs( strcat(working_folder,'\' ,cell2mat(txt1(i,20))), 'filetype', 'sfp'); 

        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,i);

        EEG = pop_saveset( EEG, 'filepath', working_folder,'filename',EEG_datasetname); 
        
        %=============================
        %ch_no=39;
        ch_no_str=cell2mat(txt2(i,21)); 
        ch_no_str=strsplit(',',ch_no_str); 
        [rc cc]=size(ch_no_str); 
        for rep_no=1:cc-1 
              ch_no=str2num(ch_no_str{rep_no});
              tmp=EEG.data(ch_no,:,:);
              RT=[];
              format_tmp=zeros(1001,EEG.trials-1); 
              for k=1:EEG.trials-1
                 format_tmp(:,k)=tmp(1,:,k);
                 RT=[RT EEG.event(2*k).latency-EEG.event(2*k-1).latency]; 
              end
              RT_sorted=sort(RT);
              mean_RT=mean(RT);middle_RT=median(RT);std_RT=std(RT);
              label_fig=strcat(cell2mat((txt1(i,22))),'_',EEG_datasetname,'_', num2str(EEG.chanlocs(ch_no).labels),'_MeanRT=',num2str(round(mean_RT)));
              figure;erpimage(format_tmp,RT_sorted,[-prestim traillength+1 Fs],label_fig,2,1,'erp','cbar','vert',0);
              savefile_erpImage=strcat(working_folder,'\ERPImage_',cell2mat((txt1(i,22))),'_',cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,2))),'_prev_stim.pdf');
              export_fig (savefile_erpImage,'-append', '-transparent')

        end 
        
        
        
        %=================================
        
    end
     
  end


end
