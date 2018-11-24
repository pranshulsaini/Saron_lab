mydata=EEG.data;
chn_no=EEG.nbchan;
Fs=1000; 
F=1:1:Fs/2; % standard Frequency Spectrum
epoch_no=EEG.trials-60;
time_bin=[0 250;250 500; 500 750; 750 1000; 1000 1250]; % start and end of each epoch bin
time_bin_no=size(time_bin);time_bin_no=time_bin_no(1); % no of bins
freq_bin=[6 8; 8 10; 10 12]; % start and end of each freq bin ( does not einclude the end point) 
freq_bin_no=size(freq_bin);freq_bin_no=freq_bin_no(1); % no of freq bins
avg_coh=zeros(1,Fs/2);
raster_mat=[];
raster_signal=[]; counter=0;



for i=1:chn_no-1
     i
     for j=i+1:chn_no
        if (i~=j)        
         for tb=1:time_bin_no          

           for k=1:epoch_no
         
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
                     raster_mat(j,i,tb,fb)=raster_mat(i,j,tb,fb);       % averaging the coherence values in frquecny bin =fb  
                     raster_mat(i,i,tb,fb)=1;
                 end
                  
         end
        end
        
     end
end



% for j=1:freq_bin_no
%   for i=1:time_bin_no
%       counter=0; 
%       for k1=1:chn_no-1
%            for k2=k1+1:chn_no
%                counter=counter+1;
%                  raster_signal(counter,i,j)= raster_mat(k1,k2,i,j);
%            end
%       end
%            
%    end
% end


 

for j=1:freq_bin_no
  for i=1:time_bin_no
       subplot(2,3,i);
       k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), 'Hz. Time frame:', num2str(time_bin(i,1)), '->' ,num2str(time_bin(i,2)), 'msec' ); 
       imagesc(raster_mat(:,:,i,j));
       title(k)
   end
     figure;
end


ras=raster_mat(:,:,1,1);
[icasig, A, W] = FASTICA (ras);
figure;imagesc(icasig);
relative_power=sum(icasig.^2')./sum(sum(icasig.^2'));
subplot(211);plot(relative_power); title(' component relative power for Raster Plot ')



thr1=0.6; thr2=0.6;
figure
k=strcat('Reconstructing the connectivty based on the original raster plot with the thersold ',num2str(thr1));
title(k);
for i=1:chn_no-1
      for j=i+1:chn_no
           
         x1=EEG.chanlocs(i).X;   y1=EEG.chanlocs(i).Y; z1=EEG.chanlocs(i).Z;label1=EEG.chanlocs(i).urchan; 
         tmp=x1; x1=-y1;y1=tmp;
         x2=EEG.chanlocs(j).X;   y2=EEG.chanlocs(j).Y; z2=EEG.chanlocs(j).Z;label2=EEG.chanlocs(j).urchan;
         tmp=x2; x2=-y2;y2=tmp;

         Wid1=  abs(ras(i,j));
         if Wid1>=thr1
            line([x1 x2],[y1 y2],[z1 z2],'Marker','.','LineStyle','-','LineWidth',Wid1*4),
            text(x1, y1, z1,num2str(label1));
            text(x2, y2,z2,num2str(label2));
           
         end

      end
end



figure;
for cmp=1:chn_no-1
icasig_selected=icasig(cmp,:);
A2=A(:,cmp);
ras_estimate=A2*icasig_selected;



  for i=1:chn_no-1
      for j=i+1:chn_no
         
         x1=EEG.chanlocs(i).X;   y1=EEG.chanlocs(i).Y; z1=EEG.chanlocs(i).Z;label1=EEG.chanlocs(i).urchan; 
         tmp=x1; x1=-y1;y1=tmp;
         x2=EEG.chanlocs(j).X;   y2=EEG.chanlocs(j).Y; z2=EEG.chanlocs(j).Z;label2=EEG.chanlocs(j).urchan;
         tmp=x2; x2=-y2;y2=tmp;
         
         subplot(4,7,cmp)
         k=strcat('comp.',num2str(cmp),' thr:' , num2str(thr2) );
         title(k); 
         Wid2=  abs(ras_estimate (i,j));
         if Wid2>=thr2 
            line([x1 x2],[y1 y2],[z1 z2],'Marker','.','LineStyle','-','LineWidth',Wid2 * 4,'color','r');
               text(x1, y1, z1,num2str(label1));
                text(x2, y2,z2,num2str(label2));
           
         end 
      end
  end
  
end 
subplot(4,7,cmp+1);plot(relative_power); title(' component relative power for Raster Plot ')

figure;
recon_comp=[1:26];
icasig_selected=icasig(recon_comp,:);
A2=A(:,recon_comp);
ras_estimate=A2*icasig_selected;
 for i=1:chn_no-1
      for j=i+1:chn_no
         
         x1=EEG.chanlocs(i).X;   y1=EEG.chanlocs(i).Y; z1=EEG.chanlocs(i).Z;label1=EEG.chanlocs(i).urchan; 
         tmp=x1; x1=-y1;y1=tmp;
         x2=EEG.chanlocs(j).X;   y2=EEG.chanlocs(j).Y; z2=EEG.chanlocs(j).Z;label2=EEG.chanlocs(j).urchan;
         tmp=x2; x2=-y2;y2=tmp;
     
         k=strcat('reconstrducted based on ',num2str(recon_comp), ' thr:' , num2str(thr2) );
         title(k); 
         Wid2=  abs(ras_estimate (i,j));
         if Wid2>=thr2 
                line([x1 x2],[y1 y2],[z1 z2],'Marker','.','LineStyle','-','LineWidth',Wid2 * 4,'color','g');
               text(x1, y1, z1,num2str(label1));
                text(x2, y2,z2,num2str(label2));
           
         end 
      end
  end
  
