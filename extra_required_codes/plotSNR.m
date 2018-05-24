% This code call SNR_CALC for estimation of SNR based on the paper from 
%Methodological Issues in Event-Related Brain Potential and Magnetic Field Studies
%Walton T. Roth, Judith M. Ford, Adolf Pfefferbaum, and Thomas R. Elbert
%here the SNR for each different experimental condition (cond) is
%calcualted and topography of SNR is depicted 
% EEGFILE_AVG: is the average EEG per condition 
% EVTFILE_AVG: is theevent file for EEGFILE_AVG
%EEGFILE: the dat file for all trials 
%EVTfile : event file for EEGFILE
% Cond: Condition of the experiment                
% Chunk: split the data into how many segments
% take attention to numbers of electrodes in .dat and .sfp. They should be
% identical
% check for the event files and sfp files for any extra line or extra
% electrodes like  FTZ....

% eegfile_avg='C:\Iman Work\tst\AVG_recon_ASDMSI_003_stim_forSOBI_ACCA_av-export.dat'
% eegfile='C:\Iman Work\tst\recon_ASDMSI_003_stim_forSOBI_ACCA.dat'
% evtfile='C:\Iman Work\tst\recon_ASDMSI_003_stim_forSOBI_ACCA.evt'
% evtfile_avg=evtfile;
% sfpfilename='C:\Iman Work\tst\forBesa_TTTTASDMSI_003_forSOBI.sfp'

% eegfile_avg='C:\Iman Work\tst\beforesobi\AVG_ASDMSI_003_forSOBI_av-export.dat'
% eegfile='C:\Iman Work\tst\beforesobi\ASDMSI_003_forSOBI.dat'
% evtfile='C:\Iman Work\tst\beforesobi\ASDMSI_003_forSOBI.evt'
% evtfile_avg=evtfile;
% sfpfilename='C:\Iman Work\tst\beforesobi\TTTTASDMSI_003_forSOBI.sfp'

% eegfile_avg='C:\Iman Work\tst\004\iman\AVG_recon_ASDMSI004_forSOBI_ACCA_av-export.dat'
% eegfile='C:\Iman Work\tst\004\iman\recon_ASDMSI004_forSOBI_ACCA.dat'
% evtfile='C:\Iman Work\tst\004\iman\recon_ASDMSI004_forSOBI_ACCA.evt'
% evtfile_avg=evtfile;
% sfpfilename='C:\Iman Work\tst\004\forBesa_ASDMSI_004_forSOBI.sfp'
% 
eegfile_avg='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_019\GOTO SOBI\SOBI_Output\ASDMSI_019_forSOBI\For SNR\recon_ASDMSI_019_forSOBI_ACCA_av-export.dat';
eegfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_019\GOTO SOBI\SOBI_Output\ASDMSI_019_forSOBI\recon_ASDMSI_019_forSOBI_ACCA.dat';
evtfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_019\GOTO SOBI\SOBI_Output\ASDMSI_019_forSOBI\recon_ASDMSI_019_forSOBI_ACCA.evt';
evtfile_avg=evtfile;
sfpfilename='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI_019\GOTO SOBI\forSNR_forBesa_ASDMSI_019_forSOBI.sfp' ;     

% eegfile_avg='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNR\ASDMSI_006_preSOBI_-200+800_Average2-export.dat';
% eegfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNR\ASDMSI_006_forSOBI-export.dat';
% evtfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNR\ASDMSI_006_forSOBI-export.evt'
% evtfile_avg=evtfile;
% sfpfilename='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNR\ASDMSI_006_forSOBI-export.sfp'      

% eegfile_avg='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNR\recon_ASDMSI_006_afterSOBI_-200+800_Average-export.dat';
% eegfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNRrecon_ASDMSI_006_forSOBI_ACCA-export.dat';
% evtfile='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNR\recon_ASDMSI_006_forSOBI_ACCA-export.evt';
% evtfile_avg=evtfile;
% sfpfilename='C:\Iman Work\CMB Projects\DATA\MSI\Subject\ASDMSI 006\GOTO SOBI\SOBI_Output\ASDMSI_006_forSOBI\for SNR\recon_ASDMSI_006_forSOBI_ACCA-export.sfp' ;     

el = readlocs(sfpfilename);
all_snr=[]; chunk=5;
for cond=11:17
    [snr_cond number_of_electrode]= SNR_CALC_CHUNKED(eegfile_avg, evtfile_avg,eegfile, evtfile,cond,chunk);
   
    all_snr=[all_snr;snr_cond];
    for i=1:chunk
        subplot(chunk,7,(cond-10)+(7*(i-1))); 
        col=(1+(number_of_electrode*(cond-11)):number_of_electrode*(cond-10));
        topoplot(all_snr(col,i), el, 'electrodes', 'off');  
    end
end

result=[];   
for i=11:17
 
    result=[result;mean(all_snr((i-10-1)*number_of_electrode+1:(i-10)*number_of_electrode))];
end
