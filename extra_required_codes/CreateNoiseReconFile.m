function CreateNoiseReconFile(xlsfilename, output)
%read the Excel file
[n,t,r] = xlsread(xlsfilename);

%creat new file with channels marked as norm changed to emg and vice versa
fid = fopen(output,'w');
if fid == -1
        fprintf(2, 'Error creating out file %s, make sure path exists.', output);
        return;
end
    for i = 1:1:size(r,1)
        if strcmp(r(i,2),'norm') == 1
            correction = 'emg';
        else
            correction = 'norm';
        end
            fprintf(fid, '%d %s\n', r{i,1}, correction);
    end
    fclose(fid)
end