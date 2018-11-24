% Aim: to extract features from ERPS on a grand average scale (over
% channels).
%
%
% Written by Manish Saggar (mishu@cs.utexas.edu) on 3/1/2011

function extractFeatureGA(csdmul, avgmul, sfpfile)

    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');   
    el = readlocs(sfpfile,'filetype','besa');

    % read CSD data from file and demean the data
    csd_dat = read_ucd_besa_mul(csdmul,1);
    csd_dat.filename = csdmul;
       
    % read AVG data from file and demean the data
    avg_dat = read_ucd_besa_mul(avgmul,2);
    avg_dat.filename = avgmul;
    
    t_0 = 200;
    t_end = 800;

    % do baseline correction for -200 to 0 msec
    t_bl_1 = t_0 - 199;
    t_bl_2 = t_0 + 0;
    
    avg_dat.db50 = avg_dat.db50 - mean(avg_dat.db50(:,t_bl_1:t_bl_2),2)*ones(1,800);
    avg_dat.db60 = avg_dat.db60 - mean(avg_dat.db60(:,t_bl_1:t_bl_2),2)*ones(1,800);
    avg_dat.db70 = avg_dat.db70 - mean(avg_dat.db70(:,t_bl_1:t_bl_2),2)*ones(1,800);
    avg_dat.db80 = avg_dat.db80 - mean(avg_dat.db80(:,t_bl_1:t_bl_2),2)*ones(1,800);

    csd_dat.db50 = csd_dat.db50 - mean(csd_dat.db50(:,t_bl_1:t_bl_2),2)*ones(1,800);
    csd_dat.db60 = csd_dat.db60 - mean(csd_dat.db60(:,t_bl_1:t_bl_2),2)*ones(1,800);
    csd_dat.db70 = csd_dat.db70 - mean(csd_dat.db70(:,t_bl_1:t_bl_2),2)*ones(1,800);
    csd_dat.db80 = csd_dat.db80 - mean(csd_dat.db80(:,t_bl_1:t_bl_2),2)*ones(1,800);
    
    % channels that we will be looking at - both hemispheres minus midline
    ch_h1 = getChs('left');
    ch_h2 = getChs('right');
    ch_all = getChs('all');

    ch_h1_order = match_str(avg_dat.labels, ch_h1);
    ch_h2_order = match_str(avg_dat.labels, ch_h2);

    avg_over_h1_db50 = mean(abs(avg_dat.db50(ch_h1_order,:)));
    avg_over_h1_db60 = mean(abs(avg_dat.db60(ch_h1_order,:)));
    avg_over_h1_db70 = mean(abs(avg_dat.db70(ch_h1_order,:)));
    avg_over_h1_db80 = mean(abs(avg_dat.db80(ch_h1_order,:)));
    
    avg_over_h2_db50 = mean(abs(avg_dat.db50(ch_h2_order,:)));
    avg_over_h2_db60 = mean(abs(avg_dat.db60(ch_h2_order,:)));
    avg_over_h2_db70 = mean(abs(avg_dat.db70(ch_h2_order,:)));
    avg_over_h2_db80 = mean(abs(avg_dat.db80(ch_h2_order,:)));
    
    csd_over_h1_db50 = mean(abs(csd_dat.db50(ch_h1_order,:)));
    csd_over_h1_db60 = mean(abs(csd_dat.db60(ch_h1_order,:)));
    csd_over_h1_db70 = mean(abs(csd_dat.db70(ch_h1_order,:)));
    csd_over_h1_db80 = mean(abs(csd_dat.db80(ch_h1_order,:)));
    
    csd_over_h2_db50 = mean(abs(csd_dat.db50(ch_h2_order,:)));
    csd_over_h2_db60 = mean(abs(csd_dat.db60(ch_h2_order,:)));
    csd_over_h2_db70 = mean(abs(csd_dat.db70(ch_h2_order,:)));
    csd_over_h2_db80 = mean(abs(csd_dat.db80(ch_h2_order,:)));    
    
    figure; set(gcf,'Position',[1 70 1920 1014]);
    subplot(2,1,1), plot([-199:600],avg_over_h1_db50,'r.-'); hold on; 
    subplot(2,1,1), plot([-199:600],avg_over_h1_db60,'g.-');
    subplot(2,1,1), plot([-199:600],avg_over_h1_db70,'b.-');
    subplot(2,1,1), plot([-199:600],avg_over_h1_db80,'k.-');
    axis([-200 600 0 3]); legend('50db','60db','70db','80db');
    title('Left Hemisphere');
    
    subplot(2,1,2), plot([-199:600],avg_over_h2_db50,'r.-'); hold on; 
    subplot(2,1,2), plot([-199:600],avg_over_h2_db60,'g.-');
    subplot(2,1,2), plot([-199:600],avg_over_h2_db70,'b.-');
    subplot(2,1,2), plot([-199:600],avg_over_h2_db80,'k.-');
    axis([-200 600 0 3]); legend('50db','60db','70db','80db');
    title('Right Hemisphere');
    suptitle(avgmul);
    
    figure; set(gcf,'Position',[1 70 1920 1014]);
    subplot(2,1,1), plot([-199:600],csd_over_h1_db50,'r.-'); hold on; 
    subplot(2,1,1), plot([-199:600],csd_over_h1_db60,'g.-');
    subplot(2,1,1), plot([-199:600],csd_over_h1_db70,'b.-');
    subplot(2,1,1), plot([-199:600],csd_over_h1_db80,'k.-');
    axis([-200 600 0 0.3]); legend('50db','60db','70db','80db');
    title('Left Hemisphere');

    
    subplot(2,1,2), plot([-199:600],csd_over_h2_db50,'r.-'); hold on; 
    subplot(2,1,2), plot([-199:600],csd_over_h2_db60,'g.-');
    subplot(2,1,2), plot([-199:600],csd_over_h2_db70,'b.-');
    subplot(2,1,2), plot([-199:600],csd_over_h2_db80,'k.-');
    axis([-200 600 0 0.3]); legend('50db','60db','70db','80db');
    title('Right Hemisphere');
    suptitle(csdmul);

end

function chs =  getChs(side)
    if strcmpi(side,'left')
        chs = {'s179','s180','s181','s182','s110','s111','s128','s135','s136','s130','s117','s103','s220','s212','s207','s202','s197',...
               'E49','s177','E33','s178','E18','s112','E6','s129','E5','s116','E14','s221','E27','s213','E43','s203','E56',...
               's174','s175','s176','s90','s99','s113','s114','s222','s115','s102','s98','s217','s214','s208','s204','s198',...
               'E48','s77','E32','s91','E17','s100','E16','s101','E15','s97','E28','s85','E44','s72','E57','s65','s73','s78','s86','s223','s93','s94','s95',...
               's96','E29','s84','s76','s71','s64','E60','s66','E47','s79','E31','s87','E30','s88','s83','E45','s70','E58','s62','s67','s74',...
               's80','s81','s82','s75','s69','s63','E59','s68','E46'};
    elseif strcmpi(side,'right')
        chs = {'s185','s186','s187','s188','s109','s126','s133','s139','s138','s131','s119','s104','s218','E26','s210','s206','s224','s196',...
            'E35','s189','E20','s190','E8','s125','E2','s132','E3','s120','E12','s173','s215','s209','E41','s200','E54','s191','s192','s193','s194','s108','s124','s123','s122','s121','E11','s105','s172','E25','s161','s205','s199','s195','E36','s154','E21','s166','E9','s107',...
            'E10','s106','s170','s171','E24','s165','s160','E40','s149','E53','s153','s148','s141','s142','s150','s155','s162','s167','E22','s168','s169','E23','s164','s159','E39','s147','E52','E50','s143','E37','s92','s151','s156','s163','s157','s158','s152','s146','s225','E38','s145',...
            's144','E50','s140','E51'};
    elseif strcmpi(side,'midline')
        chs = {'E34','s183','E19','s184','E7','s127','E1','s134','REF','s137','E4','s118','E13','s219','s216','s211','E42','s201','E55'};
    elseif strcmpi(side,'all')
        chs = cat(2,getChs('left'), getChs('right'), getChs('midline'));
    end
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