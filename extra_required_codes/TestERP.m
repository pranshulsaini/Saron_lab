%script to plot the ERP one segment at a time to see why it's so ugly%%
%Written by Chivon Powers 11/2013

function y = TestERP(eegfile, evtfile,sfpfile)

   eegfile = strtrim(eegfile);        
    
    % Variables
    hdr = readBESAsb_header(eegfile);
    
   % Create toi from Event file.
    evts = textread(evtfile,'%s');
    tmu = []; 
    trig = [];
     
    for i = 5:1:size(evts,1)
        if mod(i,5) == 0
            tmu(end+1) = str2num(cell2mat(evts(i)));
        elseif mod(i,5) == 2
            trig(end+1) = str2num(cell2mat(evts(i)));
        end
    end 
    
    % first look at all conditions together. 
    toi = [];
    trigAll = trig;
    cAll = 1;
    tmu = tmu./1e6; % converting to milliseconds.
    for i = 1:1:size(trig,2)
        toi(cAll,1) = round(tmu(i)*2048) - 3072; %Set the epoch range       
        toi(cAll,2) = round(tmu(i)*2048) + 3072;
        cAll = cAll + 1;
    end
    el = readlocs(sfpfile); %read in the number of electrodes
    sT = zeros(size(el,2), size(toi,1)*(toi(1,2)- toi(1,1)));
    blockSize = 6144;
    for t = 1:1:size(toi,1)
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2));
        tmp = x;
        % demeaning the sources from -200: 0 msec (shifted by 1 sample).
        tmp = tmp - mean(tmp(:,2662:3072),2)*ones(1,size(tmp,2));
        sT(:, (t-1)*blockSize+1:t*blockSize) = tmp;
    end
% create plot of ERP based on raw data

% first calculate the overall ERP 
    Xorig = zeros(size(el,2),blockSize);
    for t = 1:1:size(toi,1)
      %  Xorig = Xorig + A*St(:,toi(t,1):toi(t,2));  
        Xorig = Xorig + sT(:,(t-1)*blockSize+1:t*blockSize);
    end
    Xorig = Xorig./size(toi,1);
    XorigERP = mean(abs(Xorig));
    
    subplot(1,3,1), plot(XorigERP);
   % axis([0 6144 .05 2]);
    title(['All ', num2str(size(toi,1)), ' Trials']);
     
    
    %Next calculate the ERP for each condition
    blockSize = 6144;
    no_event=2;
    sT_new = zeros(size(el,2), blockSize,no_event);
    sT = zeros(size(el,2), blockSize*no_event);
    trigger_counter=zeros(no_event);

    tmu = [];
    trig = [];
    for i = 5:1:size(evts,1)
        if mod(i,5) == 0
            tmu(end+1) = str2num(cell2mat(evts(i)));
        elseif mod(i,5) == 2
            trig(end+1) = str2num(cell2mat(evts(i)));
        end
    end 
    
    for t = 1:1:size(toi,1)
        trigger_no=trig(t)/2;
        % the data should be read in nChans x Time format.
        x = readBESAsb_data(eegfile, hdr, toi(t,1)+1, toi(t,2));
        tmp = x;
    %    plot(x); pause;
        % demeaning the sources from -200: 0 msec (this is shifted by 1 sample also).
        tmp = tmp - mean(tmp(:,2662:3072),2)*ones(1,size(tmp,2));
        sT_new(:,:,trigger_no) = sT_new(:,:,trigger_no) + tmp; 
        trigger_counter(trigger_no)=trigger_counter(trigger_no)+1;
    end
    %sT = sT./size(toi,1);
    for i=1:no_event
        sT_new(:,:,i) = sT_new(:,:,i)./trigger_counter(i);
        sT=[sT sT_new(:,:,i)];
    end
    %plot the ERP for each condition
    XTrig1ERP = mean(abs(sT_new(:,:,1)));
    XTrig2ERP = mean(abs(sT_new(:,:,2)));
    subplot(1,3,2),plot(XTrig1ERP);
    axis([0 6144 .05 2]);
    title([num2str(trigger_counter(1)),' NonTarget Trials']);
    subplot(1,3,3),plot(XTrig2ERP);
    axis([0 6144 .05 2]);
    title([num2str(trigger_counter(2)),' Target Trials']);

end
