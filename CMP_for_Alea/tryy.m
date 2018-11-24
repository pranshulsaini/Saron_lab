clear all;
clc;

fid = fopen('C:\Alea\CMP_SOBI_data\CMP13011_forSOBI.evt','r');
MyText = textscan(fid,'%d %d %d','headerlines',1);
fclose(fid);

origtmu = MyText{1,1}(:); % includes artifacts' timestamp
origtrig = MyText{1,3}(:);  % includes arifacts
tmu = [];
trig = [];

start = 1;
while(origtrig(start) ~= 13)   % 13 is the code for staring of the introduction of compassion meditation
    start = start +1;
end


toiAll = [origtmu(start), origtmu(end)];
j = 0;
artifact_interval = origtmu(start);

for i = start: size(origtmu,1)

    if (origtrig(i) == 21 && origtrig(i+1) == 22)
        toiAll(end,2) = origtmu(i);
        toiAll(end+1,:) = [origtmu(i+1),origtmu(end)];
        artifact_interval = artifact_interval + (origtmu(i+1) - origtmu(i));      
    elseif (origtrig(i) == 21 && origtrig(i+1) ~= 22)
        for j = i:size(origtmu,1)
            if origtrig(j) == 22          % this j holds the index for artifact end
                break;
            end
        end
        artifact_interval_old = artifact_interval;
        artifact_interval = artifact_interval + (origtmu(j) - origtmu(i)); 
        for k = i+1:j-1
            if any(origtrig(k) == [10,13,14,15,16,17,18,19])  % for start conditions
                trig(end+1) = origtrig(k);   
                tmu(end+1) = origtmu(j) - artifact_interval; % they will now start at the end of the artifact
            elseif any(origtrig(k) == [20,23,24,25,26,27,28,29])  % for end condition
                trig(end+1) = origtrig(k);   
                tmu(end+1) = origtmu(i) - artifact_interval_old; % they will now start at the beginning of the artifact
            else
                fprintf('Invalid trigger')
            end
        end
        
        toiAll(end,2) = origtmu(i);
        toiAll(end+1,:) = [origtmu(j),origtmu(end)];
    
    elseif (i>j) && (origtrig(i) ~= 22)  % we don't want the trigger 22 and triggers between arifacts to get considered again
        trig(end+1) = origtrig(i);
        tmu(end+1) = origtmu(i) - artifact_interval + 100; % 100 is added so that it does not start from absolute 0
    end
end

tmu = tmu';
trig = trig';
        
%%
toi_diff = toiAll(:,2) - toiAll(:,1);

toi_stitch = zeros(size(toiAll));
toi_stitch(1,1) = 1;
toi_stitch(1,2) = toi_stitch(1,1) + toi_diff(1);
for i = 2:size(toiAll,1)
    toi_stitch(i,1) = toi_stitch(i-1,2) + 1;
    toi_stitch(i,2) = toi_stitch(i,1) + toi_diff(i);
end


        