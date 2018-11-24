% This code is written for creating HTML files from out files for reviewing
% previous decisions.
%
% input: xls file (saved as XLS 95 in MS Excel, otherwise Matlab can't read
% it). This is a two column file: (1) column is the path to the out file
% (full path), and (2) column is the corresponding path to the SMART
% folder. 
% NOTE: don't use headers in this file, the code assumes the first
% row is data not header.
%
%
% Written by Manish Saggar (mishu@cs.utexas.edu) on 2/3/2011

function creatHTMLfromOut(xlsfile)
    [n,t,r] = xlsread(xlsfile);
    for i = 1:1:size(r,1)
        % reading outfile
        sCat = readOutFile(r{i,1});
        
        % creating HTML
        printHTMLOldDecision(sCat, r{i,2});
    end   
end
function artifacts = readOutFile(outfile)
    data = importdata(outfile);
    artifacts = {};
    k = 1;
    if isfield(data,'textdata')
        for i = 1:1:size(data.textdata,1)
            ln = strsplit(' ', data.textdata{i});
            if strcmp(ln(2), 'emg') == 1 
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'emg'; 
                k = k + 1;
            elseif strcmp(ln(2), 'peaks') == 1 
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'peaks'; 
                k = k + 1;
            elseif strcmp(ln(2), 'eog') == 1 
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'eog'; 
                k = k + 1;
            elseif strcmp(ln(2), 'Reog') == 1
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'Reog';                 
                k = k + 1;
            elseif strcmp(ln(2), 'norm') == 1
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'norm';                 
                k = k + 1;                
            end
        end
    else
        for i = 1:1:size(data,1)
            ln = strsplit(' ', data{i});
            if strcmp(ln(2), 'emg') == 1 
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'emg'; 
                k = k + 1;
            elseif strcmp(ln(2), 'peaks') == 1 
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'peaks'; 
                k = k + 1;
            elseif strcmp(ln(2), 'eog') == 1 
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'eog'; 
                k = k + 1;
            elseif strcmp(ln(2), 'Reog') == 1
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'Reog';                 
                k = k + 1;
            elseif strcmp(ln(2), 'norm') == 1
                artifacts{k,1} = str2num(cell2mat(ln(1)));
                artifacts{k,2} = 'norm';                 
                k = k + 1;                                
            end
        end
    end    


end
function printHTMLOldDecision(sCat, pathToSmart)
    %read htmlfile
    f = fopen(strcat(pathToSmart,'/ArtifactsHTML.html'),'r');
    fnew = fopen(strcat(pathToSmart,'/ArtifactsHTMLOldDecision.html'),'w');
    while(1)
       ln = fgets(f);
       if ln == -1
           break;
       end
       if length(ln)<=1
           
           continue;
       end
           
       if strcmpi(ln(2:3),'TR')
           %remove previous checks
           lnNew = strrep(ln, 'checked','');
           %get the source number
           sN = regexp(lnNew, 'Source = \d*','match');
           sN = strsplit('=',sN{1});
           sN = str2num(strtrim(sN{2}));
           lnNew = strrep(lnNew, strcat('"',sCat{sN,2},'"'), strcat('"',sCat{sN,2},'" checked'));
           fprintf(fnew, '%s\n',lnNew);           
       else
           fprintf(fnew, '%s\n', ln);
       end
    end
    fclose(f);
    fclose(fnew);
   

end