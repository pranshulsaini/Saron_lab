% Aim: to find differences among neighboring channels.
%
%
% Written by Manish Saggar (mishu@cs.utexas.edu) on 5/27/2011

function findBridgeChannels(csdmul, avgmul, sfpfile, thresh)

    subj = strsplit('/',csdmul); subj = strrep(subj{end},'_','-');   
    el = readlocs(sfpfile,'filetype','besa');

    % read CSD data from file and demean the data
    csd_dat = read_ucd_besa_mul(csdmul,1);
    csd_dat.filename = csdmul;
       
    % read AVG data from file and demean the data
    avg_dat = read_ucd_besa_mul(avgmul,2);
    avg_dat.filename = avgmul;
    
    % get rid of fiducials
    if findstr(el(1).labels, 'fid')
        el(1:3) = [];
    end
       
    corrMat = zeros(length(el), length(el));
    
    diffBet = zeros(length(el), 8); % assuming maximum of 8 neighbors
    susp = {};
    k = 1;
    for ch = 1:1:length(el)
        nbs = getNeighbors(el(ch).labels);
        for n = 1:1:length(nbs)
            nn = match_str({el.labels}, nbs{n});
            if isempty(nn)
                continue;
            end
            diffBet(ch,n) = corr(csd_dat.db50(ch,:)', csd_dat.db50(nn,:)'); 
            if diffBet(ch,n) >= thresh
               if corrMat(ch,nn) == -1
                   continue;
               end
               susp{k} = {el(ch).labels, nbs{n},diffBet(ch,n)};
               k = k + 1;               
               corrMat(ch,nn) = -1;
               corrMat(nn,ch) = -1;
            end
        end
    end
    clc;
    fprintf(2, 'Using CSD data, found %d suspicious bridging cases (beware of duplicates), using %f as threshold for correlations. They are:\n',length(susp), thresh);
    for s = 1:1:length(susp)
        fprintf(1, '%s and its neighbor %s with corr = %f \n',cell2mat(susp{s}(1)), cell2mat(susp{s}(2)),cell2mat(susp{s}(3)));
        if s < length(susp)
            fprintf(1, 'and \n');
        end
    end

    corrMat = zeros(length(el), length(el));
    diffBet = zeros(length(el), 8); % assuming maximum of 8 neighbors
    susp = {};
    k = 1;
    for ch = 1:1:length(el)
        nbs = getNeighbors(el(ch).labels);
        for n = 1:1:length(nbs)
            nn = match_str({el.labels}, nbs{n});
            if isempty(nn)
                continue;
            end
            diffBet(ch,n) = corr(avg_dat.db50(ch,:)', avg_dat.db50(nn,:)'); 
            if diffBet(ch,n) > thresh
               if corrMat(ch,nn) == -1
                   continue;
               end 
               susp{k} = {el(ch).labels, nbs{n},diffBet(ch,n)};
               k = k + 1;
               corrMat(ch,nn) = -1;
               corrMat(nn,ch) = -1;
            end
        end
    end

    fprintf(2, 'Using AVG data, found %d suspicious bridging cases (beware of duplicates), using %f as threshold for correlations. They are:\n',length(susp), thresh);
    for s = 1:1:length(susp)
        fprintf(1, '%s and its neighbor %s with corr = %f \n',cell2mat(susp{s}(1)), cell2mat(susp{s}(2)),cell2mat(susp{s}(3)));
        if s < length(susp)
            fprintf(1, 'and \n');
        end
    end
    
    close all;
    figure; set(gcf, 'Position',[19         -21         817        1105]);
    np = length(susp);
    for s = 1:1:length(susp)
        subplot(np, 2, 2*(s-1)+1), plot(avg_dat.db50(match_str({el.labels},susp{s}(1)),:),'k--'); hold on;
        plot(avg_dat.db50(match_str({el.labels},susp{s}(2)),:),'r-');
        title(strcat('Channel(Black)-',cell2mat(susp{s}(1)),'-and-Channel(Red)-',cell2mat(susp{s}(2)),', Corr=',num2str(cell2mat(susp{s}(3)))));
        
        subplot(np, 2, 2*(s-1)+2), plot(avg_dat.db50(match_str({el.labels},susp{s}(1)),:)-...
            avg_dat.db50(match_str({el.labels},susp{s}(2)),:),'k--'); axis([1 length(avg_dat.timeline) -0.1 0.1]);
        title(strcat('Difference between Channel-',cell2mat(susp{s}(1)),'-and-',cell2mat(susp{s}(2))));
    end
    
    suptitle(strcat('Finding bridges in:',subj,'-Using threshold = ',num2str(thresh)));
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf, strrep(subj,'.mul','_bridgeCh.tiff'), 'tiff');
    close all; 
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
function neighbors = getNeighbors(ch)
    neighbors_name = {   
                    {'E2','REF','E6','E18','E7','E8'}
                    {'E8','E9','E10','E3','REF','E1'}
                    {'E2','E10','E11','E12','E4','REF'}
                    {'REF','E3','E12','E13','E14','E5'}
                    {'E6','REF','E4','E14','E15','E16'}
                    {'E18','E1','REF','E5','E16','E17'}
                    {'E19','E20','E8','E1','E18','E33'}
                    {'E20','E21','E9','E2','E1','E7','E19'}
                    {'E21','E22','E10','E2','E8'}
                    {'E9','E22','E23','E11','E3','E2'}
                    {'E10','E23','E24','25','E12','E3'}
                    {'E3','E11','E25','E26','E13','E4','E24'}
                    {'E4','E12','E26','E27','E14'}
                    {'E5','E4','E13','E27','E28','E15','E29'}
                    {'E16','E5','E14','E28','E29','E30'}
                    {'E17','E6','E5','E15','E29','E30','E31'}
                    {'E32','E18','E6','E16','E31','E30'}
                    {'E33','E7','E1','E6','E17','E32'}
                    {'E34','E35','E20','E7','E33','E49','E8','E18'}
                    {'E35','E36','E21','E8','E7','E19','E34'}
                    {'E36','E37','E22','E9','E8','E20'}
                    {'E21','E37','E38','E23','E10','E9'}
                    {'E22','E38','E39','E24','E11','E10'}
                    {'E23','E39','E40','E25','E11'}
                    {'E11','E24','E40','E41','E26','E12'}
                    {'E12','E25','E41','E42','E27','E13'}
                    {'E14','E13','E26','E42','E43','E28'}
                    {'E15','E14','E27','E43','E44','E29'}
                    {'E30','E16','E15','E28','E44','E45'}
                    {'E31','E16','E15','E29','E45','E46','E17'}
                    {'E32','E17','E16','E30','E46','E47'}
                    {'E48','E33','E18','E17','E31','E47'}
                    {'E49','E34','E19','E7','E18','E32','E48'}
                    {'E35','E20','E19','E33','E49'}
                    {'E21','E36','E21','E20','E19','E34'}
                    {'E50','E37','E21','E20','E35'}
                    {'E50','E51','E38','E22','E21','E36'}
                    {'E51','E37','E22','E23','E39'}
                    {'E38','E23','E24','E40','E52','E53'}
                    {'E39','E24','E25','E41','E53','E52','E54'}
                    {'E25','E26','E42','E54','E40','E55'}
                    {'E26','E27','E43','E56','E55','E54','E41'}
                    {'E27','E28','E44','E56','E42','E57'}
                    {'E29','E28','E43','E57','E58','E45'}
                    {'E44','E29','E30','E46','E58','E57'}
                    {'E30','E31','E47','E59','E58','E45'}
                    {'E31','E32','E48','E60','E59','E46'}
                    {'E32','E47','E60','E49'}
                    {'E34','E19','E33','E32','E48'}
                    {'E36','E37','E51'}
                    {'E38','E37','E50'}
                    {'E39','E53','E40'} 
                    {'E52','E39','E40','E41','E54'}
                    {'E41','E40','E53','E55','E42'}
                    {'E42','E41','E54','E56','E43'}
                    {'E57','E44','E43','E42','E55'}
                    {'E45','E44','E58','E56'}
                    {'E44','E46','E57'}
                    {'E60','E47','E46'}
                    {'E48','E47','E59'}
                    {'E1','E2','E3','E4','E5','E6'}
               };
           
            ch_names = {
                    'E1'
                    'E2'
                    'E3'
                    'E4'
                    'E5'
                    'E6'
                    'E7'
                    'E8'
                    'E9'
                    'E10'
                    'E11'
                    'E12'
                    'E13'
                    'E14'
                    'E15'
                    'E16'
                    'E17'
                    'E18'
                    'E19'
                    'E20'
                    'E21'
                    'E22'
                    'E23'
                    'E24'
                    'E25'
                    'E26'
                    'E27'
                    'E28'
                    'E29'
                    'E30'
                    'E31'
                    'E32'
                    'E33'
                    'E34'
                    'E35'
                    'E36'
                    'E37'
                    'E38'
                    'E39'
                    'E40'
                    'E41'
                    'E42'
                    'E43'
                    'E44'
                    'E45'
                    'E46'
                    'E47'
                    'E48'
                    'E49'
                    'E50'
                    'E51'
                    'E52'
                    'E53'
                    'E54'
                    'E55'
                    'E56'
                    'E57'
                    'E58'
                    'E59'
                    'E60'
                    'REF'
                };
         
         if isempty(match_str(ch_names', ch))
             neighbors = {};
         else
             neighbors = neighbors_name{match_str(ch_names', ch)};    
         end
end
