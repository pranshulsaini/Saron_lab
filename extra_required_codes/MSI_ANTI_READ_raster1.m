function MSI_ANTI_READ_raster1(xlsfile) 
% updtaed 12/21/2012
%function [ALLEEG,EEG,CURRENTSET]= MSI_ANTI_EEGLAB(xlsfile,percent_data,freq_bin) 
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
        cd (working_folder)
        d=dir('*.mat'); [num_file,c]=size(d);
        elec_num=(cell2mat(txt2(i,15)));
       
        %===================================================
        ALLRASTERS = zeros(elec_num, elec_num,9); 
        ALL_mean_thr=0; % only for one freq bin
        ALL_std_thr=0;
        filelist=[];
        accu_raster_mat=zeros(elec_num, elec_num,9); 
        for fn=1:num_file
           clear raster_mat mean_thr std_thr thr1 thr2
           filename=d(fn).name;
           load (filename)
           tmp_ras=raster_mat>thr1;
           accu_raster_mat=accu_raster_mat+tmp_ras;
           ALLRASTERS=ALLRASTERS+ raster_mat;
           ALL_mean_thr=ALL_mean_thr+mean_thr;
           ALL_std_thr=ALL_std_thr+std_thr;
           filelist=[filelist;{filename}];
        end
        accu_raster_mat=accu_raster_mat/num_file;
        ALLRASTERS=ALLRASTERS/num_file;
        ALL_mean_thr=ALL_mean_thr/num_file;
        ALL_std_thr=ALL_std_thr/num_file;
        comment='EVOKED'; 
       [EEG]=pop_loadset('ANTI_template_dataset.set', working_folder);
        datasetname=strcat(cell2mat((txt1(i,1))),'_',cell2mat((txt1(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
       [status]=coherence_all_raster1(EEG, datasetname, working_folder,ALLRASTERS,accu_raster_mat,ALL_mean_thr,ALL_std_thr ,prestim,filelist,comment);
        if status==1
           display(strcat('Average EVOKED Coherence Maps have been generated for datasets....',datasetname));
        end 
        
        %===============INDUCED AVERAGE================
        
        ALLRASTERS = zeros(elec_num, elec_num,9); 
        ALL_mean_thr=0; % only for one freq bin
        ALL_std_thr=0;
        filelist=[];
        accu_raster_mat=zeros(elec_num, elec_num,9); 
        for fn=1:num_file
           clear raster_mat_induced mean_thr_induced std_thr_induced thr1_induced thr2_induced
           filename=d(fn).name;
           load (filename)
           tmp_ras=raster_mat_induced>thr1_induced;
           accu_raster_mat=accu_raster_mat+tmp_ras;
           ALLRASTERS=ALLRASTERS+ raster_mat_induced;
           ALL_mean_thr=ALL_mean_thr+mean_thr_induced;
           ALL_std_thr=ALL_std_thr+std_thr_induced;
           filelist=[filelist;{filename}];
        end
        accu_raster_mat=accu_raster_mat/num_file;
        ALLRASTERS=ALLRASTERS/num_file;
        ALL_mean_thr=ALL_mean_thr/num_file;
        ALL_std_thr=ALL_std_thr/num_file;
       comment='INDUCED'; 
       [EEG]=pop_loadset('ANTI_template_dataset.set', working_folder);
        datasetname=strcat(cell2mat((txt1(i,1))),'_',cell2mat((txt1(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
       [status]=coherence_all_raster1(EEG, datasetname, working_folder,ALLRASTERS,accu_raster_mat,ALL_mean_thr,ALL_std_thr ,prestim,filelist,comment);
        if status==1
           display(strcat('Average INDUCED Coherence Maps have been generated for datasets....',datasetname));
        end 
        
    end
  end
  
end
