freq_bin=[];
mydata=EEG.data;
chn_no=EEG.nbchan;
Fs=1000; 
F=1:1:Fs/2; % standard Frequency Spectrum
epoch_no=EEG.trials-71;
time_bin=[100 200;200 300; 300 400; 400 500; 500 600; 600 700; 700 800; 800 900; 900 1000]; % start and end of each epoch bin
time_bin_no=size(time_bin);time_bin_no=time_bin_no(1); % no of bins
clc
freq_bin=input('Enter  Freq Bin (i.e. [6 8; 8 10; 10 12; 6 12; 35 45] ) you want to continue proceed with that:> ');

if isempty(freq_bin) 
   freq_bin=[6 8; 8 10; 10 12; 6 12; 35 45]; % start and end of each freq bin ( does not einclude the end point) 
end

freq_bin_no=size(freq_bin);freq_bin_no=freq_bin_no(1); % no of freq bins
avg_coh=zeros(1,Fs/2);
raster_mat=[];
raster_signal=[]; counter=0;


clc
for i=1:chn_no-1
     k=strcat('Calculating Coherence for Channel....', num2str(i));
     disp(k)
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
                 ind_chk=find(Cxy_new>1); 
                 Cxy_new(ind_chk)=1;
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

datasetname=strrep(EEG.filename,'.set','');
savefile=strcat(eeglab_evtfile,'\',datasetname,'.pdf');

% thereshold for coherence based on Ding paper
mean_thr=[];std_thr=0;
for j=1:freq_bin_no
    figure;
    scrsz = get(0,'ScreenSize'); figure('Position',[1 1 2000 2000]);
    m=0;s=0;
   for i=1:time_bin_no
       subplot(3,3,i);
       k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), 'Hz. TF:', num2str(time_bin(i,1)), '->' ,num2str(time_bin(i,2)), 'msec' ); 
       imagesc(raster_mat(:,:,i,j));
       title(k)
       m=m+mean(mean(raster_mat(:,:,i,j)));
       s=s+mean(std(raster_mat(:,:,i,j)));
   end;
   mean_thr(j)=m/time_bin_no;std_thr(j)=s/time_bin_no;
   export_fig (savefile,'-append')

     
end


%==================================================

% Drawing Coherence map


% desired_freq_bin=input('Enter Freq Bin Number (ie. 1 or 2 ,...) you want to continue proceed with that:> ');
high_node_weight=zeros(time_bin_no,freq_bin_no,chn_no);
low_node_weight=zeros(time_bin_no,freq_bin_no,chn_no);
for j=1:freq_bin_no
    figure;
    scrsz = get(0,'ScreenSize'); figure('Position',[1 1 2000 2000]);

    for i=1:time_bin_no
        thr1=mean_thr(j)+(2*std_thr(j)); % thereshold for coherence based on Ding paper
        thr2=mean_thr(j)-(2*std_thr(j)); % for finding desync regions 
        
       % k=strcat('the connectivity plot based on the original raster plot using thershold value :>  ',num2str(thr1));
        %title(k);
        ras=raster_mat(:,:,i,j);
        subplot(3,3,i);
        k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), ' TF:', num2str(time_bin(i,1)), '->' ,num2str(time_bin(i,2)), 'ms thr:  ',num2str(thr1) ); 
        title(k);
        [x_elec,y_elec] = topoplot2([],EEG.chanlocs,'electrodes','numbers');z_elec = ones(size(x_elec))*2.1;
        hold on; plot3(y_elec,x_elec,z_elec,'.r');
        for ii=1:chn_no-1
                
               for jj=i+1:chn_no
                  if ii~=jj
                    Wid1=  abs(ras(ii,jj));
                    if Wid1>=thr1
%           
                         line([y_elec(ii) y_elec(jj)],[x_elec(ii) x_elec(jj)],[z_elec(ii) z_elec(jj)],'Marker','.','LineStyle','-','LineWidth',Wid1*4),
                         high_node_weight(i,j,ii)=  high_node_weight(i,j,ii)+1;
                         high_node_weight(i,j,jj)=  high_node_weight(i,j,jj)+1;

                       
                    end
                    
                    if Wid1<=thr2
%           
                         line([y_elec(ii) y_elec(jj)],[x_elec(ii) x_elec(jj)],[z_elec(ii) z_elec(jj)],'Color','r','Marker','.','LineStyle','-','LineWidth',Wid1*4),
                         low_node_weight(i,j,ii)=  low_node_weight(i,j,ii)+1; %tb fb elecnum
                         low_node_weight(i,j,jj)=  low_node_weight(i,j,jj)+1;
                    end
                  end 
                    
                    

                end
        end

    end
    export_fig (savefile,'-append')
end
%========================================================

figure;
for fb=1:freq_bin_no
    strength_mat1=[];
    strength_mat2=[];

    for tb=1:time_bin_no
        tmp1=high_node_weight(tb,fb,:);tmp1=tmp1(:)';
        tmp2=low_node_weight(tb,fb,:);tmp2=tmp2(:)';
        
        strength_mat1(tb,:)=tmp1;
        strength_mat2(tb,:)=tmp2;
    end
       subplot(freq_bin_no,2,2*fb-1)
       imagesc(strength_mat1');xlabel('Time Bins'); ylabel('Channel No.');k=strcat('Freq: ', num2str(freq_bin(fb,1)), '->', num2str(freq_bin(j,2)), ' Hi. Coh.'); title(k);
       
       subplot(freq_bin_no,2,2*fb); 
       imagesc(strength_mat2');xlabel('Time Bins'); ylabel('Channel No.');k=strcat('Freq: ', num2str(freq_bin(fb,1)), '->', num2str(freq_bin(j,2)), ' Lo. Coh.'); title(k);
       
end



