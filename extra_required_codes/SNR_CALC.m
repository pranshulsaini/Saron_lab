
 
function y = SNR_CALC(eegfile_avg, evtfile_avg,eegfile, evtfile,cond)
%Methodological Issues in Event-Related Brain Potential and Magnetic Field Studies
%Walton T. Roth, Judith M. Ford, Adolf Pfefferbaum, and Thomas R. Elbert
%here the SNR for each different experimental condition (cond) is
%calcualted and topography of SNR is depicted 
% EEGFILE_AVG: is the average EEG per condition 
% EVTFILE_AVG: is theevent file for EEGFILE_AVG
%EEGFILE: the dat file for all trials 
%EVTfile : event file for EEGFILE
% Cond: Condition of the experiment 
%===========================================================
    eegfile = strtrim(eegfile);        
    
    % Variables
    hdr = readBESAsb_header(eegfile);
 
    
   % Creating toi from Event file.
    evts = textread(evtfile,'%s');
    tmu = []; 
    trig = [];
     
%     for i = 5:1:size(evts,1)
%         if mod(i,4) == 1
%             tmu(end+1) = str2num(cell2mat(evts(i)));
%         elseif mod(i,4) == 3
%             trig(end+1) = str2num(cell2mat(evts(i)));
%         end
%     end 
    
    for i = 5:1:size(evts,1)
        if mod(i,5) == 0
            tmu(end+1) = str2num(cell2mat(evts(i)));
        elseif mod(i,5) == 2
            trig(end+1) = str2num(cell2mat(evts(i)));
        end
    end
    % keeping all conditions together. 
    toiAll = [];
    trigAll = trig;
    cAll = 1;
    tmu = tmu./1000; % converting to milliseconds.
    for i = 1:1:size(trig,2)
        toiAll(cAll,1) = tmu(i) - 199; %from ANTI stim-locked, data shifted by 1 sample compared to APP. Check if it applies here.        
        toiAll(cAll,2) = tmu(i) + 800;
        cAll = cAll + 1;
    end
 %========================================================================   
eegfile_avg = strtrim(eegfile_avg);        
    
    % Variables
    hdr_avg = readBESAsb_header(eegfile_avg);
 
    
   % Creating toi from Event file.
    evts_avg = textread(evtfile_avg,'%s');
    tmu_avg = []; 
    trig_avg = [];
     
%     for i = 5:1:size(evts_avg,1)
%         if mod(i,4) == 1
%             tmu_avg(end+1) = str2num(cell2mat(evts_avg(i)));
%         elseif mod(i,4) == 3
%             trig_avg(end+1) = str2num(cell2mat(evts_avg(i)));
%         end
%     end 
%     
%     % keeping all conditions together. 
% %     toiAll_avg = [];
%     trigAll_avg = trig;
%     cAll_avg = 1;
%     tmu_avg = tmu_avg./1000; % converting to milliseconds.
%     for i = 1:1:size(trig_avg,2)
%         toiAll_avg(cAll_avg,1) = tmu_avg(i) - 199; %from ANTI stim-locked, data shifted by 1 sample compared to APP. Check if it applies here.        
%         toiAll_Avg(cAll_avg,2) = tmu_avg(i) + 800;
%         cAll_avg = cAll_avg + 1;
%     end
    %====================================================================================       
    x_avg = readBESAsb_data(eegfile_avg, hdr_avg, toiAll(cond-10,1), toiAll(cond-10,2));
    x_avg_max=max(x_avg');
    
    
    %x_avg=mean(x);
  
    I=find(trig==cond);
    x_cond=[];
    for i=1:size(I,2)
        x = readBESAsb_data(eegfile, hdr, toiAll(I(i),1), toiAll(I(i),2));
        x_cond=[x_cond, x];
    end 
    x_cond_max=max(x_cond');
    
    
    number_of_electrode=size(x_avg,1);
    T=size(x_avg,2);
    y=[]; 
    for ch=1:number_of_electrode
       Xt=x_avg(ch,:)/x_avg_max(ch);
       X_est=0;
       x_cond_ch=x_cond(ch,:)/x_cond_max(ch);
       for i=1:T
          X_est=X_est+Xt(i)^2;
       end
       X_est=X_est/T;
       sigma_n=0;
       for trial_no=1:size(I,2)
           tmp=x_cond_ch((trial_no-1)*1000+1:trial_no*1000);

          for t=1:T
              tmp=tmp/max(tmp);
              sigma_n=sigma_n+(Xt(t)-tmp(t))^2;
          
          end
       end
       sigma_n=sigma_n/(T*(size(I,2)-1));
       sigma_s=X_est-sigma_n/size(I,2);
       snr=sigma_s/sigma_n;   
       y=[y; snr]; 
       
       
    end
   
        
end
