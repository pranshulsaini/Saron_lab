function CreateGoodTrialsEventfileRIT(nameoffile)

UIDfilename = strcat(nameoffile, '_UniqueIDs.evt');
%eegfile = strcat(nameoffile, '_epoched_forSOBI');
inloc = '\\DSS02721-cmb-d\g\RIT\PreSOBI\';

%Read in original events 
UIDevents = textread(strcat(inloc, UIDfilename),'%s');
UID_tmu = [];
UID_trig = [];
UID_Block = [];
UID_ISI = [];
UID_RTmu = [];
OrigTrigs = [];
UID_TrialCount = [];

P1GoodTrials_trig = [];
P1GoodTrials_tmu = [];
P1GoodTrials_Block = [];
P1GoodTrials_ISI = [];
P1GoodTrials_RTmu = [];
P1GoodTrials_OGtrigs = [];
P1GoodTrials_Count = [];

P2GoodTrials_trig = [];
P2GoodTrials_tmu = [];
P2GoodTrials_Block = [];
P2GoodTrials_ISI = [];
P2GoodTrials_RTmu = [];
P2GoodTrials_OGtrigs = [];
P2GoodTrials_Count = [];

%Read in Complete Trial  Information
for j = 9:9:size(UIDevents)
    
    UID_trig(end+1) = str2num(cell2mat(UIDevents(j+2)));
    UID_tmu(end+1) = str2num(cell2mat(UIDevents(j)));
    UID_Block(end+1) = str2num(cell2mat(UIDevents(j+5)));
    UID_ISI(end+1) = str2num(cell2mat(UIDevents(j+6)));
    UID_RTmu(end+1) = str2num(cell2mat(UIDevents(j+7)));
    OrigTrigs(end+1) = str2num(cell2mat(UIDevents(j+4)));
    UID_TrialCount(end+1) = str2num(cell2mat(UIDevents(j+8)));
    
end
%Read in 1st pass Trials
Pass1filename = strcat(nameoffile, '_1stPass.evt');

outloc = '\\DSS02721-cmb-d\g\RIT\PreSOBI\';
P1events = textread(strcat(inloc, Pass1filename),'%s');
P1_trig = [];
P1_tmu = [];

for j = 6:3:size(P1events)
    P1_trig(end+1) = str2num(cell2mat(P1events(j)));
    P1_tmu(end+1) = str2num(cell2mat(P1events(j-2)));
end

for j = 1:1:size(P1_trig,2)
    
    for i = 1:1:size(UID_trig,2)
        if strcmp(num2str(UID_trig(i)),num2str(P1_trig(j))) && strcmp(num2str(UID_tmu(i)),num2str(P1_tmu(j)))
            P1GoodTrials_trig(end+1) = UID_trig(i);
            P1GoodTrials_tmu(end+1) = UID_tmu(i);
            P1GoodTrials_Block(end+1) = UID_Block(i);
            P1GoodTrials_ISI(end+1) = UID_ISI(i);
            P1GoodTrials_RTmu(end+1) = UID_RTmu(i);
            P1GoodTrials_OGtrigs(end+1) = OrigTrigs(i);
            P1GoodTrials_Count(end+1) = UID_TrialCount(i);
            
            if i < size(UID_trig,2) && strcmp(num2str(UID_trig(i+1)),'1')
                P1GoodTrials_trig(end+1) = UID_trig(i+1);
                P1GoodTrials_tmu(end+1) = UID_tmu(i+1);
                P1GoodTrials_Block(end+1) = UID_Block(i+1);
                P1GoodTrials_ISI(end+1) = UID_ISI(i+1);
                P1GoodTrials_RTmu(end+1) = UID_RTmu(i+1);
                P1GoodTrials_OGtrigs(end+1) = 1;
                P1GoodTrials_Count(end+1) = UID_TrialCount(i+1);
            end
            break
        end
        
        
        
    end
end

%Read in 2nd pass Trials
Pass2filename = strcat(nameoffile, '_BlinkPass.evt');

outloc = '\\DSS02721-cmb-d\g\RIT\PreSOBI\';
P2events = textread(strcat(inloc, Pass2filename),'%s');
P2_trig = [];
P2_tmu = [];

for j = 6:3:size(P2events)
    P2_trig(end+1) = str2num(cell2mat(P2events(j)));
    P2_tmu(end+1) = str2num(cell2mat(P2events(j-2)));
end

for j = 1:1:size(P1GoodTrials_trig,2)
    
    for i = 1:1:size(P2_trig,2)
        if strcmp(num2str(P1GoodTrials_trig(j)),num2str(P2_trig(i))) && strcmp(num2str(P1GoodTrials_tmu(j)),num2str(P2_tmu(i)))
            P2GoodTrials_trig(end+1) = P1GoodTrials_trig(j);
            P2GoodTrials_tmu(end+1) = P1GoodTrials_tmu(j);
            P2GoodTrials_Block(end+1) = P1GoodTrials_Block(j);
            P2GoodTrials_ISI(end+1) = P1GoodTrials_ISI(j);
            P2GoodTrials_RTmu(end+1) = P1GoodTrials_RTmu(j);
            P2GoodTrials_OGtrigs(end+1) = P1GoodTrials_OGtrigs(j);
            P2GoodTrials_Count(end+1) = P1GoodTrials_Count(j);
            
            if j < size(P1GoodTrials_trig,2) && strcmp(num2str(P1GoodTrials_trig(j+1)),'1')
                P2GoodTrials_trig(end+1) = P1GoodTrials_trig(j+1);
                P2GoodTrials_tmu(end+1) = P1GoodTrials_tmu(j+1);
                P2GoodTrials_Block(end+1) = P1GoodTrials_Block(j+1);
                P2GoodTrials_ISI(end+1) = P1GoodTrials_ISI(j+1);
                P2GoodTrials_RTmu(end+1) = P1GoodTrials_RTmu(j+1);
                P2GoodTrials_OGtrigs(end+1) = 1;
                P2GoodTrials_Count(end+1) = P1GoodTrials_Count(j+1);
            end
            break
        end
         
    end
    
   
    
    
end

    TargetHits = 0;
    TargetMiss = 0;
    NonTargetCRs = 0;
    NonTargetFA = 0;
    
  
for i = 1:1:size(P2GoodTrials_OGtrigs,2)
    if i == size(P2GoodTrials_OGtrigs,2)
        if P2GoodTrials_OGtrigs(i) == 2
            NonTargetFA = NonTargetFA + 1;
        elseif P2GoodTrials_OGtrigs(i) == 4
            TargetHits = TargetHits +1;
        end
    elseif P2GoodTrials_OGtrigs(i)== 2
        if P2GoodTrials_OGtrigs(i+1) == 2 || P2GoodTrials_OGtrigs(i+1) == 4
            NonTargetFA = NonTargetFA + 1;
        elseif  P2GoodTrials_OGtrigs(i+1) == 1
            NonTargetCRs = NonTargetCRs + 1;
        end
    elseif P2GoodTrials_OGtrigs(i)== 4
        if P2GoodTrials_OGtrigs(i+1) == 2 || P2GoodTrials_OGtrigs(i+1) == 4
            TargetHits = TargetHits + 1;
        elseif P2GoodTrials_OGtrigs(i+1) == 1
            TargetMiss = TargetMiss + 1;
        end
    end
    
    
end
   Total =  NonTargetFA+ NonTargetCRs+ TargetMiss+ TargetHits;
%Write out eventfiles
writeGTevt(strcat([outloc, nameoffile,'_GoodTrials.evt']), P2GoodTrials_trig, P2GoodTrials_tmu,P2GoodTrials_Block,P2GoodTrials_ISI,P2GoodTrials_RTmu, P2GoodTrials_OGtrigs, P2GoodTrials_Count);
%writeSOBIevt(strcat([inloc, nameoffile,'_epoched_forSOBI.evt']), P2GoodTrials_trig, P2GoodTrials_OGtrigs);
writeConditions(strcat([outloc, nameoffile,'_ConditionCounts.txt']), NonTargetFA, NonTargetCRs, TargetMiss, TargetHits, Total);
end


function writeGTevt(evtfile, P2GoodTrials_trig, P2GoodTrials_tmu,P2GoodTrials_Block,P2GoodTrials_ISI,P2GoodTrials_RTmu, P2GoodTrials_OGtrigs, P2GoodTrials_Count)
% creating event file
fid = fopen(evtfile,'w');
if fid == -1
    fprintf(1,'Error creating event file.\n');
    return;
end
fprintf(fid,'Tmu\tCode\tTriNo\tComnt\tBlock\tISI\tRTmu\tTrialCount\n');

for t = 1:1:size(P2GoodTrials_trig,2)
    
    fprintf(fid, '%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\n', P2GoodTrials_tmu(t), 1, P2GoodTrials_trig(t), ['Trig. ',num2str(P2GoodTrials_OGtrigs(t))], num2str(P2GoodTrials_Block(t)),P2GoodTrials_ISI(t),P2GoodTrials_RTmu(t), P2GoodTrials_Count(t));
    
end
disp('GoodTrials Event file was successfully created.');
fclose(fid);
end

function writeConditions(txtfile, NonTargetFAs, NonTargetCRs, TargetMiss, TargetHits, Total)
% creating event file
fid2 = fopen(txtfile,'w');
if fid2 == -1
    fprintf(1,'Error creating event file.\n');
    return;
end
fprintf(fid2,'FalseAlarms\tCorrectRejections\tMisses\tHits\tTotal\n');


    fprintf(fid2, '%d\t%d\t%d\t%d\t%d\n', NonTargetFAs, NonTargetCRs, TargetMiss, TargetHits, Total);
    
disp('Condition file was successfully created.');
fclose(fid2);
end

function writeSOBIevt(evtfile, P2GoodTrials_trig, P2GoodTrials_OGtrigs)
% creating event file
fid = fopen(evtfile,'w');
if fid == -1
    fprintf(1,'Error creating event file.\n');
    return;
end
fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
    epoch_time = 1500000;
for t = 1:1:size(P2GoodTrials_trig,2)
if P2GoodTrials_trig(t) > 1
     fprintf(fid, '%d\t%d\t%d\t%s\n', epoch_time, 1, P2GoodTrials_OGtrigs(t), ['Trig. ',num2str(P2GoodTrials_OGtrigs(t))]);
      epoch_time = epoch_time + 3000488;
end
end
disp('SOBI Event file was successfully created.');
fclose(fid);

end
