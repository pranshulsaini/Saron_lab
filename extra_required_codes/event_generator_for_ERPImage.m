function [latency] = event_generator_for_ERPImage2(evtfile,newevtfile,evtcode,respcode)
%evtfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_021\Goto SOBI\SOBI_Output\ASDMSI_021_forSOBI\forFieldtrip_EEGLAB\recon_ASDMSI_021_forSOBI_ACCA.evt';
% Enetr the event file from ASDMSIxxx_PTAClean_GoodTrials_RT_Triggers_edited
% Note: it should have only two columns tmu and trigger code
% Output is in mirco second scale so keep it in mind for importing in EEGLAB
latency=[]; 
evts = textread(evtfile,'%s');
    tmu = 0; 
    trig = 0;
    counter=0;
    for i = 5:1:size(evts,1)
        
        if mod(i,4) == 1 
            tmu = str2num(cell2mat(evts(i)));
        elseif mod(i,4) == 3
            trig = str2num(cell2mat(evts(i)));
        end
        if (trig==evtcode)
            counter=counter+1; 
            latency(counter)=str2num(cell2mat(evts(i+1)));
            seg_time(counter)=str2num(cell2mat(evts(i+1)))-str2num(cell2mat(evts(i-1)));
            trig=0;
        end 
            
    end 

%tmu=round(tmu./1000000);  % suppoese the Fs=1000Hz and data in evt file is coming in Tmu code trigger order and tmu is in usec
[r,c]=size(latency');
latency=round(latency); 
latency=[latency' seg_time' respcode*ones(r,1)];
Response_time=latency(:,2);

fid = fopen(strcat(newevtfile,'\MSI_ANTI_RT_Cond', cond,'EEGLAB.evt'),'w');
if fid == -1
   fprintf(1,'Error creating New Event File  \n');
end
fprintf(fid,'Response_Time\n');
for i=1:size(Response_time)
            fprintf(fid, '%s\n',num2str(Response_time(i)));
end 

    fclose(fid);    


%============================================