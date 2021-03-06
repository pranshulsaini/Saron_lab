% What yoo should do before using the code:
% 1- in the before_SOBI file , delete the last column ( Comment Column ) -
% you can do it by copy and paste all the file to Excel 
% 2- % response trig. Otherwise the code gives you an error 
% 3- Check for any blank lines in the evt files (last lines) and delete them 
% Run the code


%before_sobi_evt='C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\6030_anti_stim_-500+1100_forSOBI_edited_2.evt';
%after_sobi_evt='C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\recon_6030_anti_stim_-500+1100_forSOBI_ACCA.evt';


function [status] = MSI_ANTI_EventFileCorr2(before_sobi_evt,after_sobi_evt)
% reading the evt files 

[tmu1,code1,triger1]=textread(before_sobi_evt,'%s %s %s');
[tmu2,code2,triger2,desc2]=textread(after_sobi_evt,'%s %s %s %s');
mat_new=[];
[r1,c1]=size(tmu1);
[r2,c2]=size(tmu2);

% Checking for compatibility of the evt file
if r1-1~=(r2-1)*2 
    disp('sth wrong in the forSOBI file');
    stop
end  
 

% computing the new evt file
counter=1;
for i=2:r2   % skupping the first line -headr 
    find1=0;
    while (find1==0)  % checking for any bad trigger out of the rabge of interest
        counter=counter+1;
        tmp=str2num(cell2mat(triger1(counter)));
        if ((tmp>=11) && (tmp<=17)) ||  ((tmp>=31) && (tmp<=34)) 
            start_point=str2num(cell2mat(tmu1(counter)));
            find1=1;
        end 
    end
    
    find1=0;
    while (find1==0)
        counter=counter+1;
        tmp=str2num(cell2mat(triger1(counter)));
        if (tmp==101) ||  (tmp==444) || (tmp ==555) 
            end_point=str2num(cell2mat(tmu1(counter)));
            find1=1;
        end 
    end
     
    
    diff=end_point-start_point;
    mat_new=[mat_new;tmu2(i),code2(i),triger2(i),desc2(i);num2str(str2num(cell2mat(tmu2(i)))+diff), code1(counter),triger1(counter),desc2(i)];
end
[r3,c3]=size(mat_new);


% Writing the new evt file
fid = fopen(strcat(strrep(after_sobi_evt,'.evt',''),'_','AutoTrigerFixed.evt'),'w');
if fid == -1
   fprintf(1,'Error creating New Event File  \n');
end
fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
for i=1:r3
            fprintf(fid, '%s\t%s\t%s\t%s\n',cell2mat(mat_new(i,1)), cell2mat(mat_new(i,2)),cell2mat(mat_new(i,3)),cell2mat(mat_new(i,4))) ;
end 

    fclose(fid);    
status=1;
end
