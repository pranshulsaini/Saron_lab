close all;
ras=raster_mat(:,:,1,1);
%[icasig, A, W] = FASTICA (ras);
[w,sp] = runica (ras);
icasig=w*sp*ras;
A=inv(w*sp);
relative_power=sum(icasig.^2')./sum(sum(icasig.^2'));
subplot(211);plot(relative_power); title(' component relative power for Raster Plot ')



thr1=0.9; thr2=0.9;
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
%ras_estimate=A(:,cmp)*icasig_selected;
[ras_estimate] = icaproj(ras,A,cmp);



  for i=1:chn_no-1
      for j=i+1:chn_no
         
         x1=EEG.chanlocs(i).X;   y1=EEG.chanlocs(i).Y; z1=EEG.chanlocs(i).Z;label1=EEG.chanlocs(i).urchan; 
         tmp=x1; x1=-y1;y1=tmp;
         x2=EEG.chanlocs(j).X;   y2=EEG.chanlocs(j).Y; z2=EEG.chanlocs(j).Z;label2=EEG.chanlocs(j).urchan;
         tmp=x2; x2=-y2;y2=tmp;
         
         subplot(4,7,cmp)
         k=strcat('comp.',num2str(cmp),' thr:' , num2str(thr2) );
         title(k); 
         Wid2=  abs(ras_estimate (i,j))/max(max(ras_estimate));
         if Wid2>=thr2 
            line([x1 x2],[y1 y2],[z1 z2],'Marker','.','LineStyle','-','LineWidth',Wid2 * 4,'color','r');
               text(x1, y1, z1,num2str(label1));
                text(x2, y2,z2,num2str(label2));
           
         end 
      end
  end
  
end 
subplot(4,7,cmp+1);plot(relative_power); title(' component relative power for Raster Plot ')
[amp,ind]=findpeaks(relative_power);
figure;
recon_comp=ind;
icasig_selected=icasig(recon_comp,:);
%ras_estimate=A(:,recon_comp)*icasig_selected;
[ras_estimate] = icaproj(ras,A,recon_comp);
 
for i=1:chn_no-1
      for j=i+1:chn_no
         
         x1=EEG.chanlocs(i).X;   y1=EEG.chanlocs(i).Y; z1=EEG.chanlocs(i).Z;label1=EEG.chanlocs(i).urchan; 
         tmp=x1; x1=-y1;y1=tmp;
         x2=EEG.chanlocs(j).X;   y2=EEG.chanlocs(j).Y; z2=EEG.chanlocs(j).Z;label2=EEG.chanlocs(j).urchan;
         tmp=x2; x2=-y2;y2=tmp;
     
         k=strcat('reconstrducted based on ',num2str(recon_comp), ' thr:' , num2str(thr2) );
         title(k); 
         Wid2=  abs(ras_estimate (i,j))/max(max(ras_estimate));
         if Wid2>=thr2 
                line([x1 x2],[y1 y2],[z1 z2],'Marker','.','LineStyle','-','LineWidth',Wid2 * 4,'color','g');
               text(x1, y1, z1,num2str(label1));
                text(x2, y2,z2,num2str(label2));
           
         end 
      end
  end
  
