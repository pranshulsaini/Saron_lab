function [alpha_range,avg_psd,avg_freq,peak_freq]=IAF_Calc(EEG)
mydata=EEG.data;
chn_no=EEG.nbchan;
Fs=EEG.srate;
no_trial=size(EEG.data,2)/Fs;
%no_trial=EEG.trials;


% tmp=size(signal.trial); no_trial=tmp(2);
% tmp=size(signal.label); chn_no=tmp(1)-1; % deleting +EDF data channel
% Fs=signal.fsample;
h = spectrum.welch({'Hann','periodic'}, 2048, 0);

avg_psd=0;
for j=1:chn_no
   sum_psd=0;
   for i=1:no_trial
       x=EEG.data(j,:,i);
       Hpsd=psd(h,x,'Fs',Fs);
       sum_psd=sum_psd+(Hpsd.Data);
   end
   avg_psd=(sum_psd/no_trial)+avg_psd;
end
avg_psd=avg_psd/(chn_no); 
avg_freq=Hpsd.Frequencies;

 
% Peak Detection 
 [p,l] = findpeaks(avg_psd,'SORTSTR','descend','MINPEAKHEIGHT',.5);
 npeaks = 5;
 if length(l) < 5
            npeaks = length(l);
 end
 freqPeaks = round(avg_freq(l));
 peak_freq=[];
  for i = 1:1:npeaks   % just look at first three peaks
       if freqPeaks(i) ~= 59 && freqPeaks(i) ~= 60 && freqPeaks(i) ~= 61
              peak_freq = [peak_freq freqPeaks(i)];  
        end
  end

  % IAF estimation 
  disp(' The peak frequencies:') ; peak_freq
  figure; subplot(211); plot(avg_freq(1:65),avg_psd(1:65)); title(' Power Spectrum'); xlabel('Freq'); ylabel('PSD'); axis 'tight'
  subplot(212); semilogy(avg_freq(1:65),avg_psd(1:65)); title(' Semi Log Power Spectrum'); xlabel('Freq'); ylabel(' LOG PSD'); axis 'tight'

  f1=input('Enter  f1 --begining of the frequecny range for IAF:> ');
  f2=input('Enter  f2 --end of the frequecny range for IAF:> ');
  
  ind=find(avg_freq>=f1 & avg_freq<=f2);
  num=0;den=0;
  for i=1:length(ind)
      num=avg_psd(ind(i))*avg_freq(ind(i))+num;
      den=avg_psd(ind(i))+den;
  end
  IAP=num/den;
alpha_range=[IAP-4 IAP-2 IAP IAP+2]; % change from rounding
end