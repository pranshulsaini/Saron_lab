% Aim: This function is written for creating an envelope and graded butterfly
% plots.
%
% Input: Two filenames, CSD and AVG mul filenames.
% 
% Output: It saves Butterfly and 2-D mesh plots in the directory from which
% the function was called.
%
% Important Note: This function only takes input as one Mul file and then calls
% automatically all the other mul files. For example, you give
% subj_50db.mul, it will call other intensity files on its own (60, 70, and
% 80db), this is assuming that the other mul files are in the same
% directory.
%
% Written by Manish Saggar (mishu@cs.utexas.edu) 2010
%
% Modified by Manish Saggar (mishu@cs.utexas.edu) on 12/13/2010. Added
% functionality for user defined scalling. Now use this function as
% plotButterFlyUCD(csdmul, avgmul, usr_sc), where usr_sc is the user based
% scalling parameter. It is a 1 x 2 vector, e.g. value [-1 1] and e.g. call
% will thus be - plotButterFlyUCD(csdmul, avgmul, [-1 1]). 
% Notes: 
% (1) if you don't want to use usr_sc, just send the array empty as this
% e.g. does - plotButterFlyUCD(csdmul, avgmul, [])
%
function y = plotButterFlyUCD_editRJH(csdmul,avgmul,usr_sc_csd,usr_sc_averef)
    % read CSD data from file and demean the data
    csd_dat = read_ucd_besa_mul(csdmul,1);
    csd_dat.filename = csdmul;   
    
    % to scale data across conditions within subject
    dat2{1}.data = csd_dat.db50(:,find(csd_dat.timeline==0):find(csd_dat.timeline==300));
    dat2{2}.data = csd_dat.db60(:,find(csd_dat.timeline==0):find(csd_dat.timeline==300));
    dat2{3}.data = csd_dat.db70(:,find(csd_dat.timeline==0):find(csd_dat.timeline==300));
    dat2{4}.data = csd_dat.db80(:,find(csd_dat.timeline==0):find(csd_dat.timeline==300));
    
    if isempty(usr_sc_csd)
        sc(1) = max([max(max(dat2{1}.data)),max(max(dat2{2}.data)),...
            max(max(dat2{3}.data)), max(max(dat2{4}.data))]);
        sc(2) = min([min(min(dat2{1}.data)),min(min(dat2{2}.data)),...
            min(min(dat2{3}.data)), min(min(dat2{4}.data))]);
        sc = max(abs(sc(1)), abs(sc(2)));
        % making the scale 30% more than 0-300 ms time points for better
        % clarity.
        sc1 = sc + sc * 0.3;        
    else
        sc(1) = usr_sc_csd(1);
        sc(2) = usr_sc_csd(2);
        sc = max(abs(sc(1)), abs(sc(2)));
        % making the scale 30% more than 0-300 ms time points for better
        % clarity.
        sc1 = sc + sc * 0.3;  
    end
    csd_dat.sc = sc1;
    scrsz = get(0,'ScreenSize');  
    figure('Visible','off','Position',[1921 1 1200 1527]);
    set(gcf,'PaperPosition',[0,0,8.5,11])
    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');
    subplot(4,2,1),  plot(csd_dat.timeline, csd_dat.db50','b'); axis([-200 600 -sc1 sc1]); line([0 0],[sc1 -sc1],'Color','k','LineWidth',1); %title(strcat(subj,'-50DB'));    
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc1 sc1 sc1 sc1 sc1 sc1; -sc1 -sc1 -sc1 -sc1 -sc1 -sc1],'Color','k','LineWidth',1,'LineStyle',':');
    title('CSD','FontSize',12); text(700,(sc),'50db','FontSize',12);
    text(600,(sc1)+0.4,cell2mat(strsplit('-CSD.mul',subj)),'FontSize',14);    
    subplot(4,2,3),  plot(csd_dat.timeline, csd_dat.db60','b'); axis([-200 600 -sc1 sc1]); line([0 0],[sc1 -sc1],'Color','k','LineWidth',1);%title(strcat(subj,'-60DB'));
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc1 sc1 sc1 sc1 sc1 sc1; -sc1 -sc1 -sc1 -sc1 -sc1 -sc1],'Color','k','LineWidth',1,'LineStyle',':');
    text(700,(sc1),'60db','FontSize',12);
    subplot(4,2,5),  plot(csd_dat.timeline, csd_dat.db70','b'); axis([-200 600 -sc1 sc1]); line([0 0],[sc1 -sc1],'Color','k','LineWidth',1);%title(strcat(subj,'-70DB'));
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc1 sc1 sc1 sc1 sc1 sc1; -sc1 -sc1 -sc1 -sc1 -sc1 -sc1],'Color','k','LineWidth',1,'LineStyle',':');
    text(700,(sc1),'70db','FontSize',12);
    subplot(4,2,7),  plot(csd_dat.timeline, csd_dat.db80','b'); axis([-200 600 -sc1 sc1]); line([0 0],[sc1 -sc1],'Color','k','LineWidth',1);%title(strcat(subj,'-80DB'));
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc1 sc1 sc1 sc1 sc1 sc1; -sc1 -sc1 -sc1 -sc1 -sc1 -sc1],'Color','k','LineWidth',1,'LineStyle',':');
    text(700,(sc1),'80db','FontSize',12);

    % read AVG data from file and demean the data
    avg_dat = read_ucd_besa_mul(avgmul,2);
    avg_dat.filename = avgmul;
    
    % to scale data across conditions within subject
    dat2{1}.data = avg_dat.db50(:,find(avg_dat.timeline==0):find(avg_dat.timeline==300));
    dat2{2}.data = avg_dat.db60(:,find(avg_dat.timeline==0):find(avg_dat.timeline==300));
    dat2{3}.data = avg_dat.db70(:,find(avg_dat.timeline==0):find(avg_dat.timeline==300));
    dat2{4}.data = avg_dat.db80(:,find(avg_dat.timeline==0):find(avg_dat.timeline==300));
    
    if isempty(usr_sc_averef)
        sc(1) = max([max(max(dat2{1}.data)),max(max(dat2{2}.data)),...
            max(max(dat2{3}.data)), max(max(dat2{4}.data))]);
        sc(2) = min([min(min(dat2{1}.data)),min(min(dat2{2}.data)),...
            min(min(dat2{3}.data)), min(min(dat2{4}.data))]);
        sc = max(abs(sc(1)), abs(sc(2)));

        % making the scale 30% more than 0-300 ms time points for better
        % clarity.
        sc2 = sc + sc * 0.3;
    else
        sc(1) = usr_sc_averef(1);
        sc(2) = usr_sc_averef(2);
        sc = max(abs(sc(1)), abs(sc(2)));

        % making the scale 30% more than 0-300 ms time points for better
        % clarity.
        sc2 = sc + sc * 0.3;
    end
    
    avg_dat.sc = sc2;
    subj = strsplit('/',avgmul); subj = strrep(subj{end},'_','-');
    subplot(4,2,2),  plot(avg_dat.timeline, avg_dat.db50','b'); axis([-200 600 -sc2 sc2]); line([0 0],[sc2 -sc2],'Color','k','LineWidth',1);%title(strcat(subj,'-50DB'));    
    title('AVG','FontSize',12);
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc2 sc2 sc2 sc2 sc2 sc2; -sc2 -sc2 -sc2 -sc2 -sc2 -sc2],'Color','k','LineWidth',1,'LineStyle',':');
    subplot(4,2,4),  plot(avg_dat.timeline, avg_dat.db60','b'); axis([-200 600 -sc2 sc2]); line([0 0],[sc2 -sc2],'Color','k','LineWidth',1);%title(strcat(subj,'-60DB'));
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc2 sc2 sc2 sc2 sc2 sc2; -sc2 -sc2 -sc2 -sc2 -sc2 -sc2],'Color','k','LineWidth',1,'LineStyle',':');    
    subplot(4,2,6),  plot(avg_dat.timeline, avg_dat.db70','b'); axis([-200 600 -sc2 sc2]); line([0 0],[sc2 -sc2],'Color','k','LineWidth',1);%title(strcat(subj,'-70DB'));
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc2 sc2 sc2 sc2 sc2 sc2; -sc2 -sc2 -sc2 -sc2 -sc2 -sc2],'Color','k','LineWidth',1,'LineStyle',':');
    subplot(4,2,8),  plot(avg_dat.timeline, avg_dat.db80','b'); axis([-200 600 -sc2 sc2]); line([0 0],[sc2 -sc2],'Color','k','LineWidth',1);%title(strcat(subj,'-80DB'));
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc2 sc2 sc2 sc2 sc2 sc2; -sc2 -sc2 -sc2 -sc2 -sc2 -sc2],'Color','k','LineWidth',1,'LineStyle',':'); 
 
    print('-dtiff','-r300', strcat(subj,'_ButterFlyPlots.tiff'));
    close all;
    
    
    % printing 2D mesh plots
    scrsz = get(0,'ScreenSize');  
    figure('Visible','off','Position',[1921 1 1200 1527]);
    set(gcf,'PaperPosition',[0,0,8.5,11]);    
    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');

    csd_dat.data = csd_dat.db50;
    subplot(4,2,1),plotsmesh2d(csd_dat,csd_dat.sc); %title(strcat(subj,'-50DB'));
    set(gca,'Position',[0.13 0.7673 0.336 0.18]);
    title('CSD','FontSize',12); text(725,0,'50db','FontSize',12);
    text(675,-5,cell2mat(strsplit('-CSD.mul',subj)),'FontSize',14);    


    
    csd_dat.data = csd_dat.db60;
    subplot(4,2,3),plotsmesh2d(csd_dat,csd_dat.sc); %title(strcat(subj,'-60DB'));
    set(gca,'Position',[0.13 0.5482 0.336 0.18])
    % regional lines
    text(725,0,'60db','FontSize',12);
    
    csd_dat.data = csd_dat.db70;
    subplot(4,2,5),plotsmesh2d(csd_dat,csd_dat.sc); %title(strcat(subj,'-70DB'));
    set(gca,'Position',[0.13 0.3291 0.336 0.18])
    % regional lines
    text(725,0,'70db','FontSize',12);
    
    csd_dat.data = csd_dat.db80;
    subplot(4,2,7),plotsmesh2d(csd_dat,csd_dat.sc); %title(strcat(subj,'-80DB'));
    set(gca,'Position',[0.13 0.11 0.336 0.18])
    % regional lines
    text(725,0,'80db','FontSize',12);
    
    subj = strsplit('/',avgmul); subj = strrep(subj{end},'_','-');    
    if length(avg_dat.labels) > 60
        fprintf(1,'Found channels (%d) greater than 60\n',length(avg_dat.labels));
        numCh = 60;
    else
        numCh = length(avg_dat.labels);
        fprintf(1,'Found channels (%d) leq to 60\n',length(avg_dat.labels));
    end
    avg_dat.labels = avg_dat.labels(1:numCh);
    avg_dat.Channels = numCh;
    avg_dat.data = avg_dat.db50(1:numCh,:);    
    
    subplot(4,2,2),plotsmesh2d(avg_dat,avg_dat.sc); %title(strcat(subj,'-50DB'));
    set(gca,'Position',[0.57 0.7673 0.336 0.18])
    % regional lines
    title('AVG','FontSize',12);

    
    avg_dat.data = avg_dat.db60(1:numCh,:);
    subplot(4,2,4),plotsmesh2d(avg_dat,avg_dat.sc); %title(strcat(subj,'-60DB'));
    set(gca,'Position',[0.57 0.5482 0.336 0.18])
    % regional lines
    
    avg_dat.data = avg_dat.db70(1:numCh,:);
    subplot(4,2,6),plotsmesh2d(avg_dat,avg_dat.sc); %title(strcat(subj,'-70DB'));
    set(gca,'Position',[0.57 0.3291 0.336 0.18])
    % regional lines
    
    avg_dat.data = avg_dat.db80(1:numCh,:);
    subplot(4,2,8),plotsmesh2d(avg_dat,avg_dat.sc); %title(strcat(subj,'-80DB'));
    set(gca,'Position',[0.57 0.11 0.336 0.18]);
    % regional lines
   
    print('-dtiff','-r300', strcat(subj,'_2Dmesh.tiff'));
    close all;
end

function plotsmesh2d(dat,sc)
    
    chOrder{2} = [18    17    16    15    14    13     7     6     5     1     4     8     2     3    12     9    10    11];
    chOrder{4} = [49    33    34    19    35    20    61];
    chOrder{1} = [60    59    58    57    48    47    46    45    44    32    31    30    29    28];
    chOrder{5} = [43    27    56    42    55    26    54    41];
    chOrder{3} = [50    51    52    53    36    37    38    39    40    21    22    23    24    25];
    chOrderNames={'LT','C','RT','F','O'};
    chNames = {'E1'    'E2'    'E3'    'E4'    'E5'    'E6'    'E7'    'E8'    'E9'    'E10'    'E11'    'E12'    'E13'    'E14'...
               'E15'    'E16'    'E17'    'E18'    'E19'    'E20'    'E21'    'E22'    'E23'    'E24'    'E25'    'E26'    'E27'...
               'E28'    'E29'    'E30'    'E31'    'E32'    'E33'    'E34'    'E35'    'E36'    'E37'    'E38'    'E39'    'E40'...
               'E41'    'E42'    'E43'    'E44'    'E45'    'E46'    'E47'    'E48'    'E49'    'E50'    'E51'    'E52'    'E53'...
               'E54'    'E55'    'E56'    'E57'    'E58'    'E59'    'E60'    'REF'};
    
    i = [chOrder{1:5}]';       
    if size(dat.labels{1},2) > 3 % it means it contains an AVR in the labels
        dat.labels = strrep(dat.labels,'_AVR','');
    end
    
    [s1, s2]=match_str(chNames(i),dat.labels);      
   
    mesh(dat.timeline,[1:dat.Channels],dat.data(s2,:))
    axis([-200 600 1 dat.Channels -sc sc -sc sc]);
    colorbar;
    view(2);
    xlabel('Time(in msec)');
    set(gca,'YTickLabelMode','manual');
    set(gca,'YTick',[1:length(dat.labels)]);
    set(gca,'YTickLabel',dat.labels(s2),'FontSize',6);
    set(gca,'YDir','reverse');
    set(gcf, 'PaperPositionMode', 'auto');

    % regional lines
    fl = sum(ismember(dat.labels,{chNames{chOrder{1}}}));
    sl = fl + sum(ismember(dat.labels,{chNames{chOrder{2}}}));
    tl = sl + sum(ismember(dat.labels,{chNames{chOrder{3}}}));
    fol = tl + sum(ismember(dat.labels,{chNames{chOrder{4}}}));
    
    line([-200 600],[fl fl],[sc sc],'LineWidth',1,'Color','k');
    line([-200 600],[sl sl],[sc sc],'LineWidth',1,'Color','k');
    line([-200 600],[tl tl],[sc sc],'LineWidth',1,'Color','k');
    line([-200 600],[fol fol],[sc sc],'LineWidth',1,'Color','k');
    
    text(-280,fol+4,'O','FontSize',12);
    text(-280,tl+4,'F','FontSize',12);
    text(-280,sl+4,'RT','FontSize',12);
    text(-280,fl+4,'C','FontSize',12);
    text(-280,2+4,'LT','FontSize',12);
end
function dat = read_ucd_besa_mul(filename,type)
    dat = [];
    fid = fopen(filename, 'rt');
    hdr1 = fgetl(fid);  % contains Time pts./Channels etc. information
    tmp = strsplit('TimePoints=',hdr1,'append');
    dat.TimePoints = str2num(tmp{2}(1:4));%str2num(hdr1(14:17));
    tmp = strsplit('Channels=',hdr1,'append');
    dat.Channels = str2num(tmp{2}(1:4));%str2num(hdr1(28:29));
    hdr2 = fgetl(fid);  % contains channel labels
    if type == 1
        dat.labels = strtrim(strsplit('_csd',hdr2));%strsplit('	',hdr2);
    else
        dat.labels = strtrim(strsplit('_avr',hdr2));%strsplit('	',hdr2);
    end        
    if isempty(dat.labels{end})
        dat.labels = {dat.labels{1:end-1}};
    end
    dat.data = fscanf(fid, '%g');
    dat.data = reshape(dat.data, [dat.Channels dat.TimePoints]);
    dat.timeline = -200:1:599;%dat.data(1,:);

    % de-meaning the data
    dat.db50 = dat.data(:,1:800);
    dat.db60 = dat.data(:,801:1600);    
    dat.db70 = dat.data(:,1601:2400);
    dat.db80 = dat.data(:,2401:3200);
    dat = rmfield(dat,'data');
    fclose(fid);
end