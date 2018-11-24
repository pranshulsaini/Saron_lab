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
% Modified by Ryan Hubbard (rjhubbard@ucdavis.edu).  Added functionality
% for defining separate scales for CSD and AVEREF plots.  Replaced 1x2
% vector input with single value input; this single value is used for both the
% maximum and minimum intensities.  Changed scaling from 30% to 25% for
% increased clarity.  Increased number of possible channels to accomodate
% for virtual channel montages.
%
function y = SOBIplot_edit_RL_T1(csdmul,avgmul,usr_sc_csd,usr_sc_averef)
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
        sc3 = max(abs(sc(1)), abs(sc(2)));
        % making the scale 25% more than 0-300 ms time points for better
        % clarity.
        sc1 = sc3 + sc3 * 0.25;        
    else
        sc3 = usr_sc_csd;
        % making the scale 25% more than 0-300 ms time points for better
        % clarity.
        sc1 = sc3 + sc3 * 0.25;  
    end
    csd_dat.sc = sc1;
    scrsz = get(0,'ScreenSize');  
    figure('Visible','off','Position',[1921 1 1200 1527]);
    set(gcf,'PaperPosition',[0,0,8.5,11])
    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');
    subplot(4,2,1),  plot(csd_dat.timeline, csd_dat.db50','b'); axis([-200 600 -sc1 sc1]); line([0 0],[sc1 -sc1],'Color','k','LineWidth',1); %title(strcat(subj,'-50DB'));    
    line([-100 100 200 300 400 500; -100 100 200 300 400 500],[sc1 sc1 sc1 sc1 sc1 sc1; -sc1 -sc1 -sc1 -sc1 -sc1 -sc1],'Color','k','LineWidth',1,'LineStyle',':');
    title('CSD','FontSize',12); text(700,(sc1),'50db','FontSize',12);
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
        sc4 = max(abs(sc(1)), abs(sc(2)));

        % making the scale 25% more than 0-300 ms time points for better
        % clarity.
        sc2 = sc4 + sc4 * 0.25;
    else
        sc4 = usr_sc_averef;

        % making the scale 25% more than 0-300 ms time points for better
        % clarity.
        sc2 = sc4 + sc4 * 0.25;
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
    figure('Visible','off','Position',[1921 1 1000 scrsz(4)]);
    set(gcf,'PaperPosition',[0,0,8.5,14]);    
    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');

    csd_dat.data = csd_dat.db50;
    subplot(4,2,1),plotsmesh2d(csd_dat,sc3); %title(strcat(subj,'-50DB'));
    set(gca,'Position',[0.11 0.745 0.34 0.22]);
    title('CSD','FontSize',12); text(725,0,'50db','FontSize',12);
    text(675,-5,cell2mat(strsplit('-CSD.mul',subj)),'FontSize',14);    


    
    csd_dat.data = csd_dat.db60;
    subplot(4,2,3),plotsmesh2d(csd_dat,sc3); %title(strcat(subj,'-60DB'));
    set(gca,'Position',[0.11 0.51 0.34 0.22])
    % regional lines
    text(725,0,'60db','FontSize',12);
    
    csd_dat.data = csd_dat.db70;
    subplot(4,2,5),plotsmesh2d(csd_dat,sc3); %title(strcat(subj,'-70DB'));
    set(gca,'Position',[0.11 0.275 0.34 0.22])
    % regional lines
    text(725,0,'70db','FontSize',12);
    
    csd_dat.data = csd_dat.db80;
    subplot(4,2,7),plotsmesh2d(csd_dat,sc3); %title(strcat(subj,'-80DB'));
    set(gca,'Position',[0.11 0.04 0.34 0.22])
    % regional lines
    text(725,0,'80db','FontSize',12);
    
    subj = strsplit('/',avgmul); subj = strrep(subj{end},'_','-');    
    if length(avg_dat.labels) > 60
        fprintf(1,'Found channels (%d) greater than 60\n',length(avg_dat.labels));
        numCh = length(avg_dat.labels);
    else
        numCh = length(avg_dat.labels);
        fprintf(1,'Found channels (%d) leq to 60\n',length(avg_dat.labels));
    end
    avg_dat.labels = avg_dat.labels(1:numCh);
    avg_dat.Channels = numCh;
    avg_dat.data = avg_dat.db50(1:numCh,:);    
    
    subplot(4,2,2),plotsmesh2d(avg_dat,sc4); %title(strcat(subj,'-50DB'));
    set(gca,'Position',[0.57 0.745 0.34 0.22])
    % regional lines
    title('AVG','FontSize',12);

    
    avg_dat.data = avg_dat.db60(1:numCh,:);
    subplot(4,2,4),plotsmesh2d(avg_dat,sc4); %title(strcat(subj,'-60DB'));
    set(gca,'Position',[0.57 0.51 0.34 0.22])
    % regional lines
    
    avg_dat.data = avg_dat.db70(1:numCh,:);
    subplot(4,2,6),plotsmesh2d(avg_dat,sc4); %title(strcat(subj,'-70DB'));
    set(gca,'Position',[0.57 0.275 0.34 0.22])
    % regional lines
    
    avg_dat.data = avg_dat.db80(1:numCh,:);
    subplot(4,2,8),plotsmesh2d(avg_dat,sc4); %title(strcat(subj,'-80DB'));
    set(gca,'Position',[0.57 0.04 0.34 0.22]);
    % regional lines
   
    print('-dtiff','-r300', strcat(subj,'_2Dmesh.tiff'));
    close all;
end

function plotsmesh2d(dat,sc)
    
    
    chOrder{1} = [48    77  32  65  73  78  60, 66  47  79  62  67  74  80 59,  91  17  86  223 31, 81 82   46  75  68  69  63, 100 16  93  94  87  30, 83  84  45  76  70  71  58  64, 101 15  95  96  88  29, 89    97  28  102 98  217 14  221,    44  85  208 214 43  213,    72  57  204 198 203 56,   174   175 176 49  177 33, 90  99  178 18  112,  113   114     222 115 6   129 5,    181 182 19  184 187 188,    44  85  208 214 43  213, 103 27  220 13  216 104 219 26,    202 197 201 55  224 196];
    chOrder{2} = [212 207 211 42  206,     116 130 117 4   118 131 119 120,   110 111 128 7   127 1   109 126 133,    179 180 34  183 185 186,    2   132 3   124 123 122 121,    190 8   125 194 108,    35  189 20  191 192 193, 200 54  199 149 53,    210 209 41  161 205 40,    12  218 215 105 173 25  172 165,    106 11  170 171 164 24,    107 10  168 169 163 23, 157 158 38  152 145 146 225, 9   166 167 162 22, 92  156 37  151 143 144 51,    36    154 21  142 150 155 50];
    chOrderNames={'L','R'};
    chNames = {'E1'    'E2'    'E3'    'E4'    'E5'    'E6'    'E7'    'E8'    'E9'    'E10'    'E11'    'E12'    'E13'    'E14'...
               'E15'    'E16'    'E17'    'E18'    'E19'    'E20'    'E21'    'E22'    'E23'    'E24'    'E25'    'E26'    'E27'...
               'E28'    'E29'    'E30'    'E31'    'E32'    'E33'    'E34'    'E35'    'E36'    'E37'    'E38'    'E39'    'E40'...
               'E41'    'E42'    'E43'    'E44'    'E45'    'E46'    'E47'    'E48'    'E49'    'E50'    'E51'    'E52'    'E53'...
               'E54'    'E55'    'E56'    'E57'    'E58'    'E59'    'E60'    'REF'   's62' 's63' 's64' 's65' 's66' 's67' 's68'...
               's69' 's70' 's71' 's72' 's73' 's74' 's75' 's76' 's77' 's78' 's79' 's80' 's81' 's82' 's83' 's84' 's85' 's86' 's87'...
               's88' 's89' 's90' 's91' 's92' 's93' 's94' 's95' 's96' 's97' 's98' 's99' 's100' 's101' 's102' 's103' 's104' 's105'...
               's106' 's107' 's108' 's109' 's110' 's111' 's112' 's113' 's114' 's115' 's116' 's117' 's118' 's119' 's120' 's121'...
               's122' 's123' 's124' 's125' 's126' 's127' 's128' 's129' 's130' 's131' 's132' 's133' 's134' 's135' 's136' 's137'...
               's138' 's139' 's140' 's141' 's142' 's143' 's144' 's145' 's146' 's147' 's148' 's149' 's150' 's151' 's152' 's153'...
               's154' 's155' 's156' 's157' 's158' 's159' 's160' 's161' 's162' 's163' 's164' 's165' 's166' 's167' 's168' 's169'...
               's170' 's171' 's172' 's173' 's174' 's175' 's176' 's177' 's178' 's179' 's180' 's181' 's182' 's183' 's184' 's185'...
               's186' 's187' 's188' 's189' 's190' 's191' 's192' 's193' 's194' 's195' 's196' 's197' 's198' 's199' 's200' 's201'...
               's202' 's203' 's204' 's205' 's206' 's207' 's208' 's209' 's210' 's211' 's212' 's213' 's214' 's215' 's216' 's217'...
               's218' 's219' 's220' 's221' 's222' 's223' 's224' 's225' };
    
    i = [chOrder{1:2}]';       
    if size(dat.labels{1},2) > 3 % it means it contains an AVR in the labels
        dat.labels = strrep(dat.labels,'_AVR','');
    end
    
    [s1, s2]=match_str(chNames(i),dat.labels);      
   
    mesh(dat.timeline,[1:dat.Channels],dat.data(:,:))
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

    
    line([-200 600],[fl fl],[sc sc],'LineWidth',1,'Color','k');
    

    
    
    text(-280,sl+4,'R','FontSize',12);
    text(-280,2+4,'L','FontSize',12);
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