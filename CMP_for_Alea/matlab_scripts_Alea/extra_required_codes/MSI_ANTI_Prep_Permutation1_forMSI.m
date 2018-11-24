
function MSI_ANTI_Prep_Permutation1_forMSI(xlsfile) 

[num,txt1,txt2]=xlsread(xlsfile); % txt1 contains string and txt2 contains numeric values
[folder_no, para_no]=size(txt1);
folder_no=folder_no-1;

for i=2:folder_no+1
    flag=cell2num(txt2(i,3));
    if flag==1             % flag means the number of trials in this dataset is greater than minimum acceptable threshold
        working_folder=cell2mat(txt1(i,1));
        %code=cell2mat(txt1(i,2));
        saving_folder=cell2mat(txt1(i,2));
        cd (working_folder)
        d=dir('*.set'); [num_file,c]=size(d);
        [EEG]=pop_loadset('MSI_template_dataset.set', working_folder);
        tmp_chan=EEG.chanlocs;
        for fn=1:num_file
            if strcmp(d(fn).name,'MSI_template_dataset.set')==0
            [EEG]=pop_loadset(strcat(d(fn).name), working_folder);
            EEG.chanlocs=tmp_chan;
            trl=EEG.trials;
            setname=strcat(d(fn).name);
            len=length(strcat(d(fn).name));
            code=setname(len-13:len-8);
            if code(1)=='_' 
                   code=setname(len-12:len-7);
            end
            
            code=strrep(code,'_','');
                      
%             for ep=1:trl
%                 EEG.epoch(ep).eventtype(1)={code};
%                 EEG.event(2*ep-1).type=code;
%                 EEG.urevent(ep).type=code;
%             end
            EEG = eeg_checkset(EEG, 'eventconsistency');
            [ALLEEG EEG CURRENTSET] = eeg_store([], EEG,fn);
            [EEG, indices] = pop_epoch( EEG, {code}, [-.2 0.8]);
            EEG = pop_rmbase( EEG, [-200 0], []);
            
            new_setname=strcat(setname(1:len-4),'_epoched.set');
            EEG = pop_saveset( EEG, 'filepath', saving_folder,'filename',new_setname);
          %  setname=strcat(saving_folder,EEG.filename,strcat(d(fn).name,'_epoched','.set'));
            binname =strcat(saving_folder,code,'_bins.txt');
            bin_info2EEG(strcat(saving_folder,new_setname),binname,strcat(saving_folder,'\Bined\',new_setname));
            end 
        end
    end
        
end

               
% GND1=sets2GND('gui','bsln',[-200 0]);
% GND2=sets2GND('gui','bsln',[-200 0]);
% GRP_diff=GNDs2GRP('gui','create_difs','yes');
% load GRP_diff.grp -mat
% gui_erp(GRP)
% %GRP=tmaxGRP(GRP,1,'time_wind',[-500 800]);
% GRP=clustGRP(GRP,1,'time_wind',[-500 800],'chan_hood',.61,'thresh_p',.05);

end