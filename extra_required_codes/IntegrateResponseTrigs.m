%This is a script to integrate response triggers into the event files that
%come out of BESA and SOBI. The script uses the timing from the recon
%output and the response trigger RTs taken from event files created during
%preprocessing of each epoched CPT file.

function IntegrateResponseTrigs(xlsfilename)

% read the XLS file with the names of files that need event files
% integreated with response triggers
    [n, t, r] = xlsread(xlsfilename);
     
   
   
    % start from i = 2, since the first row is titles
    for i = 2:1:size(r,1)
        nameoffile = r{i,1};
        reconevtfile = strcat('recon_',nameoffile, '_epoched_forSOBI_CPT.evt');
        resptrigevtfile = strcat(nameoffile, '_RespTrigs.evt');
        eegfile = strcat(nameoffile, '_epoched_forSOBI');
        outloc = 'H:\CPT_SOBI_Output\';
        recondataloc = strcat(outloc, eegfile,'\', [eegfile, '_Sources_CPT\'], [eegfile, '_recon\']);
        Trigdataloc = strcat(outloc, 'TrigResponseEvtFiles\');
        

reconTimes = textread(strcat(recondataloc,reconevtfile),'%s');
recontmu = [];
recontrig = [];
for i =5:4:size(reconTimes)
    recontmu(end+1) = str2num(cell2mat(reconTimes(i)));
end

for i =7:4:size(reconTimes)
    recontrig(end+1) = str2num(cell2mat(reconTimes(i)));
end

RespTrigTimes = textread(strcat(Trigdataloc,resptrigevtfile), '%s');
RTs = [];
nontarget = 0;
target = 0;
hit = 0;
falsealarm = 0;
corrrej = 0;
miss = 0;
for i = 7:5:size(RespTrigTimes)-5
    if strcmp(RespTrigTimes(i),'2')
        nontarget = nontarget + 1;
%if i+5 < size(RespTrigTimes,1)  
        if  strcmp(RespTrigTimes(i+5),'1')
            falsealarm = falsealarm +1;
        else
            corrrej = corrrej +1;
        end
%end
    end
    
    if strcmp(RespTrigTimes(i),'4')
        target = target + 1;
%      if   i+5 < size(RespTrigTimes,1)  
        if  strcmp(RespTrigTimes(i+5),'1')
            ResponseRT = cell2mat(RespTrigTimes(i+3));
            TargetTime = cell2mat(RespTrigTimes(i-2));
            RTs(end+1) =  str2double(ResponseRT) - str2double(TargetTime);
            hit = hit + 1;
        else
            miss = miss + 1;
        end
  %  end
     end
end

writeBESAevt(strcat(recondataloc,'recon_',eegfile,'_withRespTrigs.evt'), recontmu,recontrig,RTs);
    end
end

function writeBESAevt(evtfile, recontmu, recontrig, RTs)
    % creating event file specifically for Eyes Open data.
    fid = fopen(evtfile,'w');
    if fid == -1
        fprintf(1,'Error creating event file with integrated response triggers.\n');
        return;
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
targetcounter = 0;
    for t = 1:1:size(recontrig,2)
   
       fprintf(fid, '%d\t%d\t%d\t%s\n', recontmu(t)+1.5e6, 1, recontrig(t), strcat('Trial:',num2str(t)));
       
       if strcmp(num2str(recontrig(t)) ,'4')
           targetcounter = targetcounter+1;
           fprintf(fid, '%d\t%d\t%d\t%s\n', recontmu(t) +1.5e6+RTs(targetcounter), 1, 1, strcat('Trial:',num2str(t)));
       end
    end
    fclose(fid);  
end