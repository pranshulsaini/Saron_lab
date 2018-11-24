function MSI_ANTI_Event_Correction6(xlsfile,min_trail)
  clc;
  [num,txt1,txt2]=xlsread(xlsfile);
  [subjects_no, para_no]=size(txt1);
  subjects_no=subjects_no-1;
   
  for i=2:subjects_no+1
    before_sobi_evt=strcat(txt1(i,3),'\',txt1(i,6));
    after_sobi_evt=strcat(txt1(i,3),'\',txt1(i,7));
    before_sobi_evt=cell2mat(before_sobi_evt);
    after_sobi_evt=cell2mat(after_sobi_evt);
    prev_stim=cell2num(txt2(i,23));
    [status,output] = MSI_ANTI_EventFileCorr7(before_sobi_evt,after_sobi_evt);
    if status==1
      display(strcat('Event File is Correctly Generated for dataset no.',num2str(i-1)));
      range = sprintf('P%i', i);
      output=strrep(output,strcat(txt1(i,3),'\'),'');
      xlswrite(xlsfile,output, 1, range)
    end 
   
   working_folder=cell2mat(txt1(i,3));
   EEG_datasetname=strcat(cell2mat((txt1(i,1))),'_',mat2str(cell2mat(txt2(i,2))),'_',mat2str(cell2mat(txt2(i,4))),'_',mat2str(cell2mat(txt2(i,5))),'_',mat2str(cell2mat(txt2(i,15))));
   %fixed_event=  cell2mat(strcat(txt1(i,3),'\',txt1(i,16)));
   fixed_event=  cell2mat(strcat(txt1(i,3),'\',output));
   
   if prev_stim>0 
       [latency,flag]=event_generator_for_ERPImage6_prevstim(fixed_event,EEG_datasetname,working_folder,(cell2mat(txt2(i,4))),(cell2mat(txt2(i,5))),min_trail,prev_stim);  
   else
       [latency,flag]=event_generator_for_ERPImage6(fixed_event,EEG_datasetname,working_folder,(cell2mat(txt2(i,4))),(cell2mat(txt2(i,5))),min_trail,prev_stim);  
   end     

   range = sprintf('Q%i', i);
   output={strcat(EEG_datasetname,'_RT_PrevStim_',num2str(prev_stim),'.evt')};
   xlswrite(xlsfile,output, 1, range)
   
   range = sprintf('S%i', i);
   output={flag};
   xlswrite(xlsfile,output, 1, range)
   end 
  
  display('***Event Correction Step has been performed SUCCESSFULLY!***')
end
