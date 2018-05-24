% Aim: to make movies out of CSD and AVG mul files.
%
%
%
%
% Written by Manish Saggar(mishu@cs.utexas.edu) on 2/21/2011
function make3dMovies(csdmul, avgmul, sfpfile, elecON, usr_sc_csd, usr_sc_averef)
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
        sc1 = sc3 + sc3 * 0;   %edited by RJH 6/21/11: eliminated 25% change to utilize whole colorbar.     
    else
        sc3 = usr_sc_csd;
        % making the scale 25% more than 0-300 ms time points for better
        % clarity.
        sc1 = sc3 + sc3 * 0;  %edited by RJH 6/21/11: eliminated 25% change to utilize whole colorbar.
    end
    csd_dat.sc = sc1;
    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');   

    figure; set(gcf,'Position',[1 70 1920 1014]);
%     suptitle(subj);
%     el = readlocs(sfpfile, 'filetype', 'besa');
    el = readlocs(sfpfile);
    
    
    % check whether el has 225 channels or not.
    if length(el) > 200
        el = readlocs(sfpfile, 'filetype','besa');
    end
    
    % removing missing chanels. 
    ch_missing = setdiff(1:length(el),match_str(upper({el.labels}), upper(csd_dat.labels)));
    el(ch_missing) = [];
    
    % check for fiducials and remove them
    if findstr(el(1).labels, 'fid')
        el(1:3) = [];
    end
    
    spl_name = strcat(subj,'.spl');
    
    % creating spline file
    if ~exist(spl_name,'file')
        headplot('setup', el, spl_name);
    end
    
    % loading besa style color map
    load colormapUCD; %Edited by RJH 6/14/11.  Edited colormap, replaced colormap_besa with colormapUCD.

    aviobj = avifile(strcat(subj,'_CSD_3D.avi'));
    t_0 = 100;
    t_end = 700;
    t_step = 1;
    k = 1;
    for t = t_0:t_step:t_end
        subplot(2,4,1), headplot(csd_dat.db50(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'cbar',0,'colormap',colormapUCD); %edited RJH 6/23/11. Turned off lighting.
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[53.9571   46.0581  -16.4099]);
        subplot(2,4,5), headplot(csd_dat.db50(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'cbar',0,'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[10.3933  -69.1149  -33.7858]);
        text(140, 120, strcat('CSD - 50db, t=',num2str(t-200)),'FontSize',24);
        
        subplot(2,4,2), headplot(csd_dat.db60(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'colormap',colormapUCD); 
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[6.3342   -5.9256  -34.4722]);
        subplot(2,4,6), headplot(csd_dat.db60(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[ -9.2620   -3.6657  -29.6687]);
        text(140, 120, strcat('CSD - 60db, t=',num2str(t-200)),'FontSize',24);
  
        subplot(2,4,3), headplot(csd_dat.db70(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'colormap',colormapUCD); 
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-1.6192  -30.4659  -32.8277]);
        subplot(2,4,7), headplot(csd_dat.db70(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-17.7742   17.7347  -28.0192]);
        text(140, 120, strcat('CSD - 70db, t=',num2str(t-200)),'FontSize',24);
  
        subplot(2,4,4), headplot(csd_dat.db80(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'colormap',colormapUCD); 
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-8.6030  -52.8358  -34.4720]);
        subplot(2,4,8), headplot(csd_dat.db80(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-csd_dat.sc,csd_dat.sc],'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-26.3898   36.9865  -30.5462]);
        text(140, 120, strcat('CSD - 80db, t=',num2str(t-200)),'FontSize',24);
        text(-420, 360, subj,'FontSize',24);
        aviobj = addframe(aviobj, getframe(gcf));
    end
    close all;
    aviobj = close(aviobj);    
    
%     save('avi_obj_csd','aviobj');
%     
%     % converting to mpeg
%     aa = aviread(strcat(subj,'_CSD_3D.avi'));
%     mpgwrite(aa,colormap('jet'),strcat(subj,'_CSD_3D.mpeg'));
%     eval(cell2mat(strcat('rm',{' '},subj,'_CSD_3D.avi')));
%     clear aa;
%     
       
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
        sc2 = sc4 + sc4 * 0; %Edited by SA 8/22/11: eliminated 25% change to utilize whole color bar
    else
        sc4 = usr_sc_averef;

        % making the scale 25% more than 0-300 ms time points for better
        % clarity.
        sc2 = sc4 + sc4 * 0; %Edited by SA 8/22/11: eliminated 25% change to utilize whole color bar
    end
    
    avg_dat.sc = sc2;
    subj = strsplit('/',avgmul); subj = strrep(subj{end},'_','-');    
    close all;
    figure; set(gcf,'Position',[1 70 1920 1014]);    
    aviobj = avifile(strcat(subj,'_AVG_3D.avi'));
    
    for t = t_0:t_step:t_end
        subplot(2,4,1), headplot(avg_dat.db50(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'cbar',0,'colormap',colormapUCD); 
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[53.9571   46.0581  -16.4099]);
        subplot(2,4,5), headplot(avg_dat.db50(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'cbar',0,'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[10.3933  -69.1149  -33.7858]);
        text(140, 120, strcat('AVG - 50db, t=',num2str(t-200)),'FontSize',24);
        
        subplot(2,4,2), headplot(avg_dat.db60(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'colormap',colormapUCD); 
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[6.3342   -5.9256  -34.4722]);
        subplot(2,4,6), headplot(avg_dat.db60(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[ -9.2620   -3.6657  -29.6687]);
        text(140, 120, strcat('AVG - 60db, t=',num2str(t-200)),'FontSize',24);
  
        subplot(2,4,3), headplot(avg_dat.db70(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'colormap',colormapUCD); 
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-1.6192  -30.4659  -32.8277]);
        subplot(2,4,7), headplot(avg_dat.db70(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-17.7742   17.7347  -28.0192]);
        text(140, 120, strcat('AVG - 70db, t=',num2str(t-200)),'FontSize',24);
  
        subplot(2,4,4), headplot(avg_dat.db80(:,t),spl_name,'view',[70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'colormap',colormapUCD); 
        set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-8.6030  -52.8358  -34.4720]);
        subplot(2,4,8), headplot(avg_dat.db80(:,t),spl_name,'view',[-70 38],'electrodes',elecON,'maplimits', [-avg_dat.sc,avg_dat.sc],'colormap',colormapUCD); 
        set(gca,'FontSize',28); set(gca,'CameraViewAngle', [5.15]); set(gca,'CameraTarget',[-26.3898   36.9865  -30.5462]);
        text(140, 120, strcat('AVG - 80db, t=',num2str(t-200)),'FontSize',24);
        text(-420, 360, subj,'FontSize',24);
        
        aviobj = addframe(aviobj, getframe(gcf));

    end
    close all;
    aviobj = close(aviobj);    
    save('avi_obj_avg','aviobj');
    
%     % converting to mpeg
%     aa = aviread(strcat(subj,'_AVG_3D.avi'));
%     mpgwrite(aa,colormap('jet'),strcat(subj,'_AVG_3D.mpeg'));
%     eval(cell2mat(strcat('rm',{' '},subj,'_AVG_3D.avi')));
%     clear aa;    

    
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