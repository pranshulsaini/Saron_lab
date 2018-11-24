% Aim: to make movies out of CSD and AVG mul files.
%
%
%
%
% Written by Manish Saggar(mishu@cs.utexas.edu) on 2/21/2011
function make2dMovies(csdmul, avgmul, sfpfile, elecON, usr_sc_csd, usr_sc_averef)
    % read CSD data from file and demean the data
    csd_dat = read_ucd_besa_mul(csdmul,1);
    csd_dat.filename = csdmul;
    
    % check for REF'' label and fix it.
    if match_str(csd_dat.labels, 'REF''')
        csd_dat.labels(match_str(csd_dat.labels,'REF''')) = {'REF'};
    end
    
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
    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');   
       
    % read AVG data from file and demean the data
    avg_dat = read_ucd_besa_mul(avgmul,2);
    avg_dat.filename = avgmul;
    
    % check for REF'' label and fix it.
    if match_str(avg_dat.labels, 'REF''')
        avg_dat.labels(match_str(avg_dat.labels,'REF''')) = {'REF'};
    end
    
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
    
    figure; set(gcf,'Position',[1 313 1920 771]);
    el = readlocs(sfpfile); 
    % check whether el has 225 channels or not.
    if length(el) > 200
        el = readlocs(sfpfile, 'filetype','besa');
    end
    
    aviobj = avifile(strcat(subj,'CSD-n-AVG-2D-full.avi'));
    ch_missing = setdiff(1:length(el),match_str(upper({el.labels}), upper(csd_dat.labels)));
    el(ch_missing) = [];
    
    for t = 100:1:700
        subplot(2,4,5), topoplot(csd_dat.db50(:,t), el, 'electrodes', elecON, 'maplimits', [-csd_dat.sc,csd_dat.sc]);
        title(strcat('CSD - 50db t=',num2str(t-200)),'FontSize',20); colorbar;
        subplot(2,4,6), topoplot(csd_dat.db60(:,t), el, 'electrodes',elecON, 'maplimits', [-csd_dat.sc,csd_dat.sc]);
        title('CSD - 60db','FontSize',20);colorbar;
        subplot(2,4,7), topoplot(csd_dat.db70(:,t), el, 'electrodes',elecON, 'maplimits', [-csd_dat.sc,csd_dat.sc]);
        title('CSD - 70db','FontSize',20);colorbar;
        subplot(2,4,8), topoplot(csd_dat.db80(:,t), el, 'electrodes',elecON, 'maplimits', [-csd_dat.sc,csd_dat.sc]);
        title('CSD - 80db','FontSize',20);colorbar;
        
        subplot(2,4,1), topoplot(avg_dat.db50(:,t), el, 'electrodes',elecON, 'maplimits', [-avg_dat.sc,avg_dat.sc]);
        title(strcat('AVG - 50db t=',num2str(t-200)),'FontSize',20);colorbar;
        subplot(2,4,2), topoplot(avg_dat.db60(:,t), el, 'electrodes',elecON, 'maplimits', [-avg_dat.sc,avg_dat.sc]);
        title('AVG - 60db','FontSize',20);colorbar;
        subplot(2,4,3), topoplot(avg_dat.db70(:,t), el, 'electrodes',elecON, 'maplimits', [-avg_dat.sc,avg_dat.sc]);
        title('AVG - 70db','FontSize',20);colorbar;
        subplot(2,4,4), topoplot(avg_dat.db80(:,t), el, 'electrodes',elecON, 'maplimits', [-avg_dat.sc,avg_dat.sc]);
        title('AVG - 80db','FontSize',20);colorbar;
        text(-3.7,-2.25, subj,'FontSize',24);
        
        aviobj = addframe(aviobj, getframe(gcf));
        

    end
    close all;
    aviobj = close(aviobj);    
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