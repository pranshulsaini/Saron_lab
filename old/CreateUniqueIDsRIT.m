%This is a script to take the original triggers from the RIT EEG event
%files and create semi-unique identifiers that allow us to track each
%individual trial through processing and analysis. 

%For the RIT task, Trigger sequences indicate the following behavioral
%conditions: 2,1 = CorrectRejection 4,1 = Miss 2,2 or 2,4 = False Alarm 4,2
%= Hit 4,4 should not happen because the experiment is programmed to never have 2 hit
%trials back to back.



function CreateUniqueIDsRIT(nameoffile)

        OrigTrigFile = strcat(nameoffile, '_aux_NoBad_AvgRef_OrigTrigs.evt');
      
        inloc = 'D:\Stroop\Temp\Batch_step\';
        outloc = 'D:\Stroop\Temp\Batch_step\';
                
Origtmu = [];
Origtrigs = [];
OriginalEvents = textread(strcat(inloc,OrigTrigFile),'%s');
for i = 25:5:size(OriginalEvents)
    Origtmu(end+1) = str2num(cell2mat(OriginalEvents(i)));
    Origtrigs(end+1) = str2num(cell2mat(OriginalEvents(i+2)));
end
%Origtmu = OriginalEvents(25:5:end);
%Origtrigs = OriginalEvents(27:5:end);
Uniqtrigs = [];
UniqtrigResp = [];
TrigBlock = [];
NonTargetID = 100;
TargetID = 500;
TBlock = 1;
NonTBlock = 1;
ISIs = [];
RTmu = [];
TrialCount = [];
Trialcounter = 0;
for j = 1:1:size(Origtrigs,2)
    
    if strcmp(num2str(Origtrigs(j)),'2')
        if j == 1
            ISIs(end+1) = str2num(cell2mat(OriginalEvents(25)))-str2num(cell2mat(OriginalEvents(20)));

        elseif j>1 && strcmp(num2str(Origtrigs(j-1)),'1')
            ISIs(end+1) = Origtmu(j) - Origtmu(j-2);
        elseif j > 1 
            ISIs(end+1) = Origtmu(j) - Origtmu(j-1);
        end
        
        TrialCount(end+1) =  Trialcounter + 1; 
         Trialcounter = Trialcounter +1; 
        NonTargetID = NonTargetID+1;
        Uniqtrigs(end+1) = NonTargetID;
        UniqtrigResp(end+1) = NonTargetID;
        TrigBlock(end+1) = NonTBlock;
        if j<size(Origtrigs,2) && strcmp(num2str(Origtrigs(j+1)),'1')
            RTmu(end+1) = Origtmu(j+1) - Origtmu(j);
        else
            RTmu(end+1) = 0;
        end
    elseif strcmp(num2str(Origtrigs(j)),'4')
        if j == 1
            ISIs(end+1) = str2num(cell2mat(OriginalEvents(25)))-str2num(cell2mat(OriginalEvents(20)));
        elseif strcmp(num2str(Origtrigs(j-1)),'1')
            ISIs(end+1) = Origtmu(j) - Origtmu(j-2);
        else
            ISIs(end+1) = Origtmu(j) - Origtmu(j-1);
        end

        TrialCount(end+1) =  Trialcounter + 1;
         Trialcounter = Trialcounter +1; 
        TargetID = TargetID+1;
        Uniqtrigs(end+1)=TargetID;
        UniqtrigResp(end+1) = TargetID;
        TrigBlock(end+1) = TBlock;
         if j<size(Origtrigs,2) && strcmp(num2str(Origtrigs(j+1)),'1')
            RTmu(end+1) = Origtmu(j+1) - Origtmu(j);
        else
            RTmu(end+1) = 0;
        end
    elseif strcmp(num2str(Origtrigs(j)),'1')
        TrialCount(end+1) =  Trialcounter;
        UniqtrigResp(end+1) = 1;
        TrigBlock(end+1) = TrigBlock(end);
        RTmu(end+1) = 0;
        ISIs(end+1) = 0;
    end;
    
    if NonTargetID == 316
        NonTargetID = 100;
        NonTBlock = NonTBlock +1;
    end
    if TargetID == 524
        TargetID = 500; 
        TBlock = TBlock + 1;
    end
      
end
writeUniqueIDevt(strcat(['\\DSS02721-CMB-D\G\RIT\PreSOBI\',nameoffile,'_UniqueIDs.evt']), Uniqtrigs,UniqtrigResp, Origtmu,Origtrigs,TrigBlock, ISIs, RTmu, TrialCount)
end

function writeUniqueIDevt(evtfilename, Uniqtrigs,UniqtrigResp, Origtmu, Origtrigs,TrigBlock, ISIs, RTmu, TrialCount)
    %Opening file to write the new events into
    fid = fopen(evtfilename,'w');
    if fid == -1
        fprintf(1,'Error creating event file with integrated response triggers.\n');
        return;
    end
    fprintf(fid,'Tmu\tCode\tTriNo\tComnt\tBlock\tISI\tRTmu\tTrialCount\n');
    for t = 1:1:size(UniqtrigResp,2)
       fprintf(fid, '%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\n', Origtmu(t), 1, UniqtrigResp(t), ['Trig. ',num2str(Origtrigs(t))], num2str(TrigBlock(t)), ISIs(t), RTmu(t), TrialCount(t));
    end
    fclose(fid);  
end