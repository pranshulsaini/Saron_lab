function ERPIMAGE_GRAND(xlsfile,condition_list,phenotype_list,save_location) 
%%
%function ERPIMAGE_GRAND(xlsfile,condition_list,save_location) 
%condition_list=[11 12 13 14 15 16 17];
[num,txt1,txt2]=xlsread(xlsfile); % txt1 contains string and txt2 contains numeric values
[subjects_no, para_no]=size(txt1);
subjects_no=subjects_no-1; % considering the header of the xlsfile
%%

for ph_lst=1:length(phenotype_list)

 for cond_lst=1:length(condition_list)

  ALLTRAILS = eeg_emptyset(); % for concatenating all connditions for IAF estimation
  counter_set=0;
  close all; 
  for i=2:subjects_no+1
    flag=cell2num(txt2(i,19)); 
    stim_trig=cell2num(txt2(i,4));
    phenotype=cell2mat((txt1(i,22))); pheno_cmp=strcmp(phenotype, phenotype_list(ph_lst));
    if (flag==1 && stim_trig==condition_list(cond_lst) && pheno_cmp==1)          
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
        EEG_datasetname=strcat(cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))),'.set'); 
        ch_no_str=cell2mat(txt2(i,21)); 
        EEG = pop_loadset(EEG_datasetname , working_folder);
        
        % GFP Normalization
        EEG.data=rmbase(EEG.data,0,100:200);
        avg_EEG=zeros(EEG.nbchan,EEG.pnts);
        for t=1:EEG.trials
            avg_EEG=avg_EEG+EEG.data(:,:,t);
        end
        avg_EEG=avg_EEG/EEG.trials;
        [gfp,gd] = eeg_gfp(avg_EEG',1);
        
        for t=1:EEG.pnts
            EEG.data(:,t,:)=EEG.data(:,t,:)/gfp(t);
        end
               
         
        
        
        
        % de-meaning each channel by substracting the average of all trials
        % in a channel from each time point in each channel
%         avg=[];vari=[];
%         for n=1:EEG.nbchan
%             tmp_data=[];
%             for t=1:EEG.trials
%                tmp_data=[tmp_data EEG.data(n,:,t)]; 
%             end
%             avg=[avg mean(tmp_data)];
%             vari=[vari var(tmp_data)];
%             temp_EEG=(EEG.data(n,:,:)-avg(n))/vari(n);
%             EEG.data(n,:,:)=temp_EEG;
%         end

        
        
        % Saving data in the merged dataset 
        
        
        if counter_set==1
            [ALLTRAILS]=EEG;
        else 
            [ALLTRAILS]=pop_mergeset(ALLTRAILS,EEG);

        end
        
    end
  end
      
  %EEG_Merged_name=strcat('MERGED_',cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
  EEG_Merged_name=strcat('MERGED_',cell2mat(phenotype_list(ph_lst)),'_',num2str(condition_list(cond_lst)),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
  ALLTRAILS = pop_saveset( ALLTRAILS, 'filepath', save_location,'filename',EEG_Merged_name);

% ERPIMAGE 
  EEG=ALLTRAILS;
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
              label_fig=strcat('MERGED_',cell2mat(phenotype_list(ph_lst)),'_',num2str(condition_list(cond_lst)),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))),'_', num2str(EEG.chanlocs(ch_no).labels));
             % figure;erpimage(format_tmp,RT_sorted,[-prestim traillength+1 Fs],label_fig,10,1,'erp','cbar','vert',0);
             figure;erpimage(format_tmp,RT_sorted,[0 traillength+1 Fs],label_fig,10,1,'erp','cbar','vert',0);
              pdfname=strcat('MERGED_',cell2mat(phenotype_list(ph_lst)),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
              savefile_erpImage=strcat(save_location,'\ERPImage_',pdfname,'.pdf');
              export_fig (savefile_erpImage,'-append', '-transparent');

  end 
 end
ph_lst;
end
%%
        
