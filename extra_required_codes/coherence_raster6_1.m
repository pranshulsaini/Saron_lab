function [status]=coherence_raster6_1(EEG,datasetname, eeglab_evtfile, percent_data,freq_bin, vol_cond)
status=0;
elec_map=[4 12 20 1 5 9 13 17 21 25 2 6 10 14 18 22 26 3 7 11 15 19 23 27 8 16 24];
%======================================
X_loc=[];Y_loc=[];Z_loc=[]; e_label=[];
for i=1:EEG.nbchan
   X_loc=[X_loc; EEG.chanlocs(i).X];
   Y_loc=[Y_loc; EEG.chanlocs(i).Y];
   Z_loc=[Z_loc; EEG.chanlocs(i).Z];
   e_label=[e_label {EEG.chanlocs(i).labels}];
end

[e_label,dist_elec]=elec_distance(e_label,X_loc,Y_loc,Z_loc,0,0,0.0359); % (0,0,0.359) is the center of head inthe besa head model with radius 85mm
dist_elec=dist_elec*100 % coverting to cm
%=======================================


elec_map_rev=(1:27);
mydata=EEG.data;
chn_no=EEG.nbchan;
Fs=1000; 
F=1:1:Fs/2; % standard Frequency Spectrum
k=strcat(' The total number of trials is : ',num2str(EEG.trials),' . How many percentage of the trails do you want to use for the coherence analysis [0 = Zero Percent...100= 100%]> ');
disp(k);
%percent_data=input(k);
epoch_no=(EEG.trials)*percent_data/100; epoch_no=round(epoch_no);
time_bin=[100 200;200 300; 300 400; 400 500; 500 600; 600 700; 700 800; 800 900; 900 1000]; % start and end of each epoch bin
time_bin_no=size(time_bin);time_bin_no=time_bin_no(1); % no of bins
clc
%freq_bin=input('Enter  Freq Bin (i.e. [6 8; 8 10; 10 12; 6 12; 35 45] ) you want to continue proceed with that:> ');

if isempty(freq_bin) ==1
   freq_bin=[6  14; 25 40]; % start and end of each freq bin ( does not einclude the end point) 
end

freq_bin_no=size(freq_bin);freq_bin_no=freq_bin_no(1); % no of freq bins
avg_coh=zeros(1,Fs/2);
raster_mat=[];
raster_signal=[]; counter=0;
raster_mat_avg_coh=[];
 

for i=1:chn_no-1
     k=strcat('Calculating Coherence for Channel....', num2str(i));
     disp(k)
     for j=i+1:chn_no
        if (i~=j)
         for tb=1:time_bin_no          
           avg_coh=zeros(1,length(F));   
           for k=1:epoch_no
         
                 x=mydata(i,:,k);
                 y=mydata(j,:,k);
                 x_bin=x(time_bin(tb,1)+1:time_bin(tb,2));
                 y_bin=y(time_bin(tb,1)+1:time_bin(tb,2));
                 [Cxy,f] = mscohere(x_bin,y_bin,[],[],[],Fs);
                 Cxy_new=spline(f,Cxy,F); % for having the same length of coherence for all loops
                 ind_chk=find(Cxy_new>1); 
                 
                 if isempty(ind_chk)==0 
                     Cxy_new(ind_chk)=1;
                 end
                 
                 % reduced coherence based on Paul L. Nunez paper
               if vol_cond(1)==1  
                 a=vol_cond(2); 
                 reduced_coh=exp((1-dist_elec(i,j))/a);
                 Cxy_new=Cxy_new-reduced_coh;
                 if Cxy_new==0
                     Cxy_new=0;
                 end
               end
                 
                 avg_coh=avg_coh+Cxy_new; % this is the sum for coherence btw channel i and j for all epochs in the time bin (tb) and for all freq 
                        
           end
                 avg_coh=avg_coh./epoch_no; % % this is the average for coherence btw channel i and j for all epochs in the time bin (tb) and for all freq
                 for fb=1:freq_bin_no
                     ind=find(F>= freq_bin(fb,1) & F<freq_bin(fb,2)); % finding indices for each frequecny bin 
                     new_i=find(elec_map==i);new_j=find(elec_map==j);
                     mm= mean(avg_coh(ind));       % averaging the coherence values in frquecny bin =fb  
                     if mm>1
                         disp('wrong')
                     end
                     raster_mat(new_i,new_j,tb,fb)=mm;
                     raster_mat(new_j,new_i,tb,fb)=raster_mat(new_i,new_j,tb,fb);       % averaging the coherence values in frquecny bin =fb  
                     raster_mat(new_i,new_i,tb,fb)=1;
                 end
                  
         end
        end
        
     end
end
% Removing Auto-Coherence values 
for i=1:chn_no
  new_i=find(elec_map==i);new_j=find(elec_map==i);
  raster_mat(new_i,new_j,:,:)=0; % exp((1-dist_elec(chn_no,chn_no))/a);
end




%datasetname=strrep(EEG.filename,'.set','');
savefile=strcat(eeglab_evtfile,'\',datasetname,'.pdf');

elec_label_y={};
for l=chn_no:-1:1
   tmp=EEG.chanlocs(elec_map(l)).labels;
   tmp=strrep(tmp,'_RFR','');
   elec_label_y=[elec_label_y {EEG.chanlocs(elec_map(l)).labels}];
end

elec_label_x={};
for l=1:chn_no
   tmp=EEG.chanlocs(elec_map(l)).labels;
   tmp=strrep(tmp,'_RFR','');
   elec_label_x=[elec_label_x {tmp}];
end


%==========================new thresholding 
% thereshold for coherence based on Ding paper

time_bin_avg=[100 1000]; % start and end of each epoch bin
time_bin_no_avg=size(time_bin_avg);time_bin_no_avg=time_bin_no_avg(1); % no of bins

for i=1:chn_no-1
     k=strcat('Calculating Coherence Threshold....', num2str(i));
     disp(k)
     for j=i+1:chn_no
        if (i~=j)
         for tb=1:time_bin_no_avg          
           avg_coh=zeros(1,length(F));   
           for k=1:epoch_no
         
                 x=mydata(i,:,k);
                 y=mydata(j,:,k);
                 x_bin=x(time_bin_avg(tb,1)+1:time_bin_avg(tb,2));
                 y_bin=y(time_bin_avg(tb,1)+1:time_bin_avg(tb,2));
                 [Cxy,f] = mscohere(x_bin,y_bin,[],[],[],Fs);
                 Cxy_new=spline(f,Cxy,F); % for having the same length of coherence for all loops
                 ind_chk=find(Cxy_new>1); 
                 
                 if isempty(ind_chk)==0 
                     Cxy_new(ind_chk)=1;
                 end
                 
                 % reduced coherence based on Paul L. Nunez paper
               if vol_cond(1)==1  
                 a=vol_cond(2); 
                 reduced_coh=exp((1-dist_elec(i,j))/a);
                 Cxy_new=Cxy_new-reduced_coh;
                 if Cxy_new==0
                     Cxy_new=0;
                 end
               end
                 
                 avg_coh=avg_coh+Cxy_new; % this is the sum for coherence btw channel i and j for all epochs in the time bin (tb) and for all freq 
                        
           end
                 avg_coh=avg_coh./epoch_no; % % this is the average for coherence btw channel i and j for all epochs in the time bin (tb) and for all freq
                 for fb=1:freq_bin_no
                     ind=find(F>= freq_bin(fb,1) & F<freq_bin(fb,2)); % finding indices for each frequecny bin 
                     new_i=find(elec_map==i);new_j=find(elec_map==j);
                     mm= mean(avg_coh(ind));       % averaging the coherence values in frquecny bin =fb  
                     if mm>1
                         disp('wrong')
                     end
                     raster_mat_avg_coh(new_i,new_j,tb,fb)=mm;
                     raster_mat_avg_coh(new_j,new_i,tb,fb)=raster_mat(new_i,new_j,tb,fb);       % averaging the coherence values in frquecny bin =fb  
                     raster_mat_avg_coh(new_i,new_i,tb,fb)=1;
                 end
                  
         end
        end
        
     end
end
% Removing Auto-Coherence values 
for i=1:chn_no
  new_i=find(elec_map==i);new_j=find(elec_map==i);
  raster_mat_avg_coh(new_i,new_j,:,:)=0; % exp((1-dist_elec(chn_no,chn_no))/a);
end

mean_thr=[];std_thr=0;
for j=1:freq_bin_no
       
       mean_thr(j)=mean(mean(raster_mat_avg_coh(:,:,:,j)));
       std_thr(j)=mean(std(raster_mat_avg_coh(:,:,:,j)));
end       
%==============================================
% thereshold for coherence based on Ding paper
%mean_thr=[];std_thr=0;
for j=1:freq_bin_no
    figure;
    scrsz = get(0,'ScreenSize'); figure('Position',[1 1 2000 2000]);
    %m=0;s=0;
    clf
    im=imread('C:\Iman Work\CMB Projects\Yukari Project\Electrodes Map\electrode map 150 dpi flattened.tif');
    axes('position',[0,.6,.1,0.63])
    imshow(im)
    hold on;
    for i=1:time_bin_no
       
       subplot(3,3,i);
       %suptitle('Raster plot of Coherence Elevation During Experiment');
       k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), 'Hz. TF:', num2str(time_bin(i,1)), '->' ,num2str(time_bin(i,2)), 'msec' ); 
       imagesc(raster_mat(:,:,i,j));
       set(gca,'YTick',[1:length(elec_map)],'YTickLabel',elec_label_x)
       set(gca,'XTick',[1:length(elec_map)],'XTickLabel',elec_label_x)
       title(k)
     %  m=m+mean(mean(raster_mat(:,:,i,j)));
      % s=s+mean(std(raster_mat(:,:,i,j)));
   end;
   %mean_thr(j)=m/time_bin_no;std_thr(j)=s/time_bin_no;
     export_fig (savefile,'-append', '-transparent')

     
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
        suptitle('Topography of Coherence Elevation During Experiment(Blue lines: High Significant Connections. Red lines : Low Significant Connections)');
        k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), ' TF:', num2str(time_bin(i,1)), '->' ,num2str(time_bin(i,2)), 'ms thr:  ',num2str(thr1) ); 
        title(k);
        [x_elec,y_elec] = topoplot2([],EEG.chanlocs,'electrodes','numbers');z_elec = ones(size(x_elec))*2.1;
        hold on; plot3(y_elec,x_elec,z_elec,'.r');
        for ii=1:chn_no-1
                
               for jj=i+1:chn_no
                  if ii~=jj
                    new_ii=elec_map(ii);new_jj=elec_map(jj);
                    Wid1=  abs(ras(ii,jj));
                    if Wid1>=thr1
%           
                         line([y_elec(new_ii) y_elec(new_jj)],[x_elec(new_ii) x_elec(new_jj)],[z_elec(new_ii) z_elec(new_jj)],'Marker','.','LineStyle','-','LineWidth',Wid1*4),
                         high_node_weight(i,j,ii)=  high_node_weight(i,j,ii)+1;
                         high_node_weight(i,j,jj)=  high_node_weight(i,j,jj)+1;

                       
                    end
                    
                    if Wid1<=thr2
%           
                         line([y_elec(new_ii) y_elec(new_jj)],[x_elec(new_ii) x_elec(new_jj)],[z_elec(new_ii) z_elec(new_jj)],'Color','r','Marker','.','LineStyle','-','LineWidth',Wid1*4),
                         low_node_weight(i,j,ii)=  low_node_weight(i,j,ii)+1; %tb fb elecnum
                         low_node_weight(i,j,jj)=  low_node_weight(i,j,jj)+1;
                    end
                  end 
                    
                    

                end
        end

    end
     export_fig (savefile,'-append', '-transparent')
end
%========================================================

figure;
scrsz = get(0,'ScreenSize'); figure('Position',[1 1 2000 2000]);
suptitle('Number of Connections for each Electrode During Experment');
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
       imagesc(strength_mat1');xlabel('Time Bins'); ylabel('Channel Label');k=strcat('Freq: ', num2str(freq_bin(fb,1)), '->', num2str(freq_bin(j,2)), ' Hi. Coh.'); title(k);
       set(gca,'YTick',[1:length(elec_map)],'YTickLabel',elec_label_x);
       set(gca,'XTick',[1:time_bin_no],'XTickLabel',time_bin(:,1));
       

       subplot(freq_bin_no,2,2*fb); 
       imagesc(strength_mat2');xlabel('Time Bins'); ylabel('Channel Label.');k=strcat('Freq: ', num2str(freq_bin(fb,1)), '->', num2str(freq_bin(j,2)), ' Lo. Coh.'); title(k);
       set(gca,'YTick',[1:length(elec_map)],'YTickLabel',elec_label_x);
       set(gca,'XTick',[1:time_bin_no],'XTickLabel',time_bin(:,1));


end
     export_fig (savefile,'-append', '-transparent')


     
%==============================writing log and mat files
logfile=strcat(eeglab_evtfile,'\',datasetname,'.log1');
fid = fopen(logfile,'w');
star_txt='*****************************************************************************************************************************';
if fid == -1
   fprintf(1,'Error creating New Event File  \n');
end
fprintf(fid,'%s\t%s\n', 'Dataset File Name:', datasetname);
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'PDF File Name for all coherence maps:', savefile);
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Frequency Bin(s):', mat2str(freq_bin));
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Time Bin(s):', mat2str(time_bin));
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Variable Name which Contains Means of Threshold For Different Frequencies Bin(s):', 'mean_thr');
fprintf(fid,'%s\t%s\n', 'mean_thr=:', mat2str(mean_thr));
fprintf(fid,'%s\t%s\n', 'Variable Name which Contains Std of Threshold For Different Frequencies Bin(s):', 'std_thr');
fprintf(fid,'%s\t%s\n', 'std_thr=:', mat2str(std_thr));
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Note: Threshold for significant coonection (coherence) for each frequency bin is Calculated as  :', 'thr1=mean_thr+(2*std_thr) and thr2=mean_thr-(2*std_thr)');
fprintf(fid,'%s\t%s\n', 'thr1=', mat2str(thr1));
fprintf(fid,'%s\t%s\n', 'thr2=', mat2str(thr2));
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Variable Name which Contains number of lines with HIGH coherence (above thr1) connected to each electrode during the experiment time: ', 'high_node_weight');
fprintf(fid,'%s\t%s\n', 'Variable Name which Contains number of lines with high coherence (above thr1) connected to each electrode during the experiment time: ', 'low_node_weight');
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Variable Name which Contains Raster Plots information:', 'raster_mat');
fprintf(fid,'%s\t%s\n', 'raster_mat is a 4 dimensional matrix: raster_mat(i,j,tb,fb) ', 'i= elec_num j=elec_num tb=time_bin_num tf=freq_bin_num');
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Actual Numeber of Trials=', num2str(EEG.trials));
fprintf(fid,'%s\t%s\n', 'Percentage of the Trials used in this analysis=', num2str(percent_data));
fprintf(fid,'%s\n', star_txt);
fclose(fid);


matfile=strcat(eeglab_evtfile,'\',datasetname,'.mat');
save(matfile, 'freq_bin', 'time_bin','raster_mat','datasetname','mean_thr','std_thr','thr1','thr2','high_node_weight','low_node_weight');
close all;
status=1;
end