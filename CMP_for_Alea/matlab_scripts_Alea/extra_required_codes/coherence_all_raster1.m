function [status]=coherence_all_raster1(EEG, datasetname,eeglab_evtfile,raster_mat,accu_raster_mat,mean_thr,std_thr ,prestim,filelist,comment)
% updtaed 01/03/2012
status=0;
elec_map=[4 12 20 1 5 9 13 17 21 25 2 6 10 14 18 22 26 3 7 11 15 19 23 27 8 16 24]; % basedhe on the scalp segmentation and the order of appearence 
[subject_num ss]=size(filelist);
%======================================
% Calculating inter-electrode distannce
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
 title_for_fig=strrep(datasetname,'_',' '); 

elec_map_rev=(1:27);
mydata=EEG.data;
chn_no=EEG.nbchan;
Fs=1000; 
F=1:1:Fs/2; % standard Frequency Spectrum

%percent_data=input(k);
% you can modify the time resolution here
time_bin=[100 200;200 300; 300 400; 400 500; 500 600; 600 700; 700 800; 800 900; 900 1000]; % start and end of each epoch bin
time_bin_no=size(time_bin);time_bin_no=time_bin_no(1); % no of bins
clc
%freq_bin=input('Enter  Freq Bin (i.e. [6 8; 8 10; 10 12; 6 12; 35 45] ) you want to continue proceed with that:> ');


freq_bin=[6  14]; % start and end of each freq bin ( does not einclude the end point) 


freq_bin_no=size(freq_bin);freq_bin_no=freq_bin_no(1); % no of freq bins
avg_coh=zeros(1,Fs/2);
raster_signal=[]; counter=0;
raster_mat_avg_coh=[]; % for thresholding 

%datasetname=strrep(EEG.filename,'.set','');
savefile=strcat(eeglab_evtfile,'\',datasetname,'_',comment,'.pdf'); % creating PDF flile 

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


%========================== thresholding 
% thereshold for coherence based on Ding paper

%==============================================
% thereshold raster plots based on Ding paper
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
       k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), 'Hz. TF:', num2str(time_bin(i,1)-prestim), '->' ,num2str(time_bin(i,2)-prestim), 'msec, ', title_for_fig ); 
       imagesc(raster_mat(:,:,i,j));
       set(gca,'YTick',[1:length(elec_map)],'YTickLabel',elec_label_x)
       set(gca,'XTick',[1:length(elec_map)],'XTickLabel',elec_label_x)
       title(k)
       colorbar('EastOutside')
     %  m=m+mean(mean(raster_mat(:,:,i,j)));
      % s=s+mean(std(raster_mat(:,:,i,j)));
   end;
   %mean_thr(j)=m/time_bin_no;std_thr(j)=s/time_bin_no;
     export_fig (savefile,'-append', '-transparent')
  
end


% Accumulative map
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
       k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), 'Hz. TF:', num2str(time_bin(i,1)-prestim), '->' ,num2str(time_bin(i,2)-prestim), 'msec, ', title_for_fig ); 
       imagesc(accu_raster_mat(:,:,i,j));
       set(gca,'YTick',[1:length(elec_map)],'YTickLabel',elec_label_x)
       set(gca,'XTick',[1:length(elec_map)],'XTickLabel',elec_label_x)
       title(k)
       colorbar('EastOutside')
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
%         suptitle('Topography of Coherence Elevation During Experiment(Blue lines: High Significant Connections. Red lines : Low Significant Connections)');
        k=strcat('Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), ' TF:', num2str(time_bin(i,1)-prestim), '->' ,num2str(time_bin(i,2)-prestim), 'ms thr:  ',num2str(thr1),', ',title_for_fig ); 
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
                         line([y_elec(new_ii) y_elec(new_jj)],[x_elec(new_ii) x_elec(new_jj)],[z_elec(new_ii) z_elec(new_jj)],'Color','r','Marker','.','LineStyle','-','LineWidth',Wid1*4),
                         high_node_weight(i,j,ii)=  high_node_weight(i,j,ii)+1;
                         high_node_weight(i,j,jj)=  high_node_weight(i,j,jj)+1;

                       
                    end
                    
                    if Wid1<=thr2 
%           
                         line([y_elec(new_ii) y_elec(new_jj)],[x_elec(new_ii) x_elec(new_jj)],[z_elec(new_ii) z_elec(new_jj)],'Color','b','Marker','.','LineStyle','-','LineWidth',Wid1*4),
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
% drawing coherence map for accumulative plot
for j=1:freq_bin_no
    figure;
    scrsz = get(0,'ScreenSize'); figure('Position',[1 1 2000 2000]);

    for i=1:time_bin_no
        thr1=1/subject_num; %mean_thr(j)+(2*std_thr(j)); % thereshold for coherence based on Ding paper
        %thr2=mean_thr(j)-(2*std_thr(j)); % for finding desync regions 
        
       % k=strcat('the connectivity plot based on the original raster plot using thershold value :>  ',num2str(thr1));
        %title(k);
        ras=accu_raster_mat(:,:,i,j);
        
        subplot(3,3,i);
%         suptitle('Topography of Coherence Elevation During Experiment(Blue lines: High Significant Connections. Red lines : Low Significant Connections)');
        k=strcat('ACCU Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), ' TF:', num2str(time_bin(i,1)-prestim), '->' ,num2str(time_bin(i,2)-prestim), 'ms thr:  ',num2str(thr1),', ',title_for_fig ); 
        title(k);
        [x_elec,y_elec] = topoplot2([],EEG.chanlocs,'electrodes','numbers');z_elec = ones(size(x_elec))*2.1;
        hold on; plot3(y_elec,x_elec,z_elec,'.r');
        for ii=1:chn_no-1
                
               for jj=i+1:chn_no
                  if ii~=jj
                    new_ii=elec_map(ii);new_jj=elec_map(jj);
                    Wid1=  abs(ras(ii,jj));
                    if Wid1>thr1
%           
                         line([y_elec(new_ii) y_elec(new_jj)],[x_elec(new_ii) x_elec(new_jj)],[z_elec(new_ii) z_elec(new_jj)],'Color','r','Marker','.','LineStyle','-','LineWidth',Wid1*9),
                                                
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
suptitle(strcat('Number of Connections for each Electrode During Experment for ', title_for_fig));
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
       set(gca,'XTick',[1:time_bin_no],'XTickLabel',time_bin(:,1)-prestim);
       colorbar('EastOutside')

       subplot(freq_bin_no,2,2*fb); 
       imagesc(strength_mat2');xlabel('Time Bins'); ylabel('Channel Label.');k=strcat('Freq: ', num2str(freq_bin(fb,1)), '->', num2str(freq_bin(j,2)), ' Lo. Coh.'); title(k);
       set(gca,'YTick',[1:length(elec_map)],'YTickLabel',elec_label_x);
       set(gca,'XTick',[1:time_bin_no],'XTickLabel',time_bin(:,1)-prestim);
       colorbar('EastOutside')

end
     export_fig (savefile,'-append', '-transparent')


     
%==============================writing log and mat files

max_dist=max(max(dist_elec));
bin_connection=20;
dist_bin=max_dist/bin_connection; 
connection_dist_mat=zeros(bin_connection,time_bin_no,freq_bin_no);
for j=1:freq_bin_no
    figure;
    scrsz = get(0,'ScreenSize'); figure('Position',[1 1 2000 2000]);
    suptitle(strcat('Percentage of Connections per inter-electrode distance ', title_for_fig));

    for i=1:time_bin_no
        ras=accu_raster_mat(:,:,i,j);
            
        subplot(3,3,i);
%       
        connection_no=zeros(1,bin_connection);
        for ii=1:chn_no-1
               
               for jj=i+1:chn_no 
                  if ii~=jj
                   new_i=find(elec_map==ii);new_j=find(elec_map==jj); % new_i and new_j are indices of the electrodes based on their apperance in the raster plot
   
                    if ras(ii,jj)>0                    
                     if (dist_elec(new_i,new_j)>0 && dist_elec(new_i,new_j)<=dist_bin)
                        connection_no (1)=connection_no (1)+1;
                     
                     elseif (dist_elec(new_i,new_j)>dist_bin && dist_elec(new_i,new_j)<= 2*dist_bin)
                         
                         connection_no (2)=connection_no (2)+1;
                     elseif (dist_elec(new_i,new_j)>2*dist_bin && dist_elec(new_i,new_j)<= 3*dist_bin)
                         connection_no (3)=connection_no (3)+1;
                     elseif (dist_elec(new_i,new_j)>3*dist_bin && dist_elec(new_i,new_j)<= 4*dist_bin)
                         connection_no (4)=connection_no (4)+1;
                     elseif (dist_elec(new_i,new_j)>4*dist_bin && dist_elec(new_i,new_j)<= 5*dist_bin)
                         connection_no (5)=connection_no (5)+1;
         
                     elseif (dist_elec(new_i,new_j)>5*dist_bin && dist_elec(new_i,new_j)<= 6*dist_bin)
                         connection_no (6)=connection_no (6)+1;
                     elseif (dist_elec(new_i,new_j)>6*dist_bin && dist_elec(new_i,new_j)<= 7*dist_bin)
                         connection_no (7)=connection_no (7)+1;
                         
                      elseif (dist_elec(new_i,new_j)>7*dist_bin && dist_elec(new_i,new_j)<= 8*dist_bin)
                         connection_no (8)=connection_no (8)+1; 

                      elseif (dist_elec(new_i,new_j)>8*dist_bin && dist_elec(new_i,new_j)<= 9*dist_bin)
                         connection_no (9)=connection_no (9)+1; 
                          
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 10*dist_bin)
                         connection_no (10)=connection_no (10)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 11*dist_bin)
                         connection_no (11)=connection_no (11)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 12*dist_bin)
                         connection_no (12)=connection_no (12)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 13*dist_bin)
                         connection_no (13)=connection_no (13)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 14*dist_bin)
                         connection_no (14)=connection_no (14)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 15*dist_bin)
                         connection_no (15)=connection_no (15)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 16 *dist_bin)
                         connection_no (16)=connection_no (16)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 17*dist_bin)
                         connection_no (17)=connection_no (17)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 18*dist_bin)
                         connection_no (18)=connection_no (18)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 19*dist_bin)
                         connection_no (19)=connection_no (19)+1;
                      elseif (dist_elec(new_i,new_j)>9*dist_bin && dist_elec(new_i,new_j)<= 20*dist_bin)
                         connection_no (20)=connection_no (20)+1;
                         
                         
                         
                         
                         
                         
                         
                     end
                    end
                      
                  end 
                    
                     

                end
        end
      %  cs=spline(1:10,connection_no,1:.25:10); plot(cs) 
          x_bin=dist_bin:dist_bin:bin_connection*dist_bin;
          connection_no=cumsum(connection_no);
          connection_no=connection_no/max(connection_no);
          plot(x_bin,connection_no); 
          ylabel('Percentage of Connections ');
          xlabel('Inter Electrode Distance (cm) ')
          k=strcat('Coonections Freq: ', num2str(freq_bin(j,1)), '->', num2str(freq_bin(j,2)), ' TF:', num2str(time_bin(i,1)-prestim), '->' ,num2str(time_bin(i,2)-prestim), ', ' ); 
          title(k);
          
    connection_dist_mat(:,i,j)=connection_no;
    end
     export_fig (savefile,'-append', '-transparent')
end



logfile=strcat(eeglab_evtfile,'\',datasetname,'.log1');
fid = fopen(logfile,'w');
star_txt='*****************************************************************************************************************************';
if fid == -1
   fprintf(1,'Error creating New Event File  \n');
end
fprintf(fid,'%s\t%s\n', 'Results for Condition:', datasetname);
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'PDF File Name for all coherence maps:', savefile);
fprintf(fid,'%s\n', star_txt);
fprintf(fid,'%s\t%s\n', 'Frequency Bin(s): The frequecny you selected in the preliminary step ; otherwise ', mat2str(freq_bin));
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

[subject_num ss]=size(filelist);
fprintf(fid,'%s\t%s\n', 'number of subjects included in this study :', mat2str(subject_num));
fprintf(fid,'%s\n', 'filelist :');
for i=1:subject_num
    txt=filelist(i); txt=cell2mat(txt);txt=mat2str(txt);
    fprintf(fid,'%s\n',  txt);
    
end


fprintf(fid,'%s\n', star_txt);
fclose(fid);


matfile=strcat(eeglab_evtfile,'\',datasetname,'');
save(matfile, 'freq_bin', 'time_bin','raster_mat','datasetname','mean_thr','std_thr','thr1','thr2','high_node_weight','low_node_weight','accu_raster_mat','connection_dist_mat');
close all;
status=1;
end