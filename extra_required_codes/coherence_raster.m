mydata=EEG.data; 
Fs=1000; 
F=1:1:Fs/2; % standard Frequency Spectrum
epoch_no=EEG.trials;
time_bin=[0 200;200 400; 400 600; 600 800; 800 1000]; % start and end of each epoch bin
time_bin_no=size(time_bin);time_bin_no=time_bin_no(1); % no of bins
freq_bin=[4 7; 8 12; 18 30]; % start and end of each freq bin ( does not einclude the end point) 
freq_bin_no=size(freq_bin);freq_bin_no=freq_bin_no(1); % no of freq bins
avg_coh=zeros(1,Fs/2);
raster_mat=[];



for i=1:EEG.nbchan-100
     i
     for j=1:EEG.nbchan-100
        for tb=1:time_bin_no          

           for k=1:epoch_no-60
         
                 x=mydata(i,:,k);
                 y=mydata(j,:,k);
                 x_bin=x(time_bin(tb,1)+1:time_bin(tb,2));
                 y_bin=y(time_bin(tb,1)+1:time_bin(tb,2));
                 [Cxy,f] = mscohere(x_bin,y_bin,[],[],[],Fs);
                 Cxy_new=spline(f,Cxy,F); % for having the same length of coherence for all loops
                 avg_coh=avg_coh+Cxy_new; % this is the sum for coherence btw channel i and j for all epochs in the time bin (tb) and for all freq 
                        
           end
                 avg_coh=avg_coh./epoch_no; % % this is the average for coherence btw channel i and j for all epochs in the time bin (tb) and for all freq
                 for fb=1:freq_bin_no
                     ind=find(F>= freq_bin(fb,1) & F<freq_bin(fb,2)); % finding indices for each frequecny bin 
                     raster_mat(i,j,tb,fb)= mean(avg_coh(ind));       % averaging the coherence values in frquecny bin =fb  
                                  
                 end
                  
        end
        
    end
end

for i=1:time_bin_no
    for j=1:freq_bin_no
       subplot(2,2,j);
       k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), 'Hz. Time frame:', num2str(time_bin(i,1)), '->' ,num2str(time_bin(i,2)), 'msec' ); 
       imagesc(raster_mat(:,:,i,j));
       title(k)

    end
     pause;
end

       