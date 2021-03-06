% What yoo should do before using the code:
% 1- Delete the first line (titles/headers) in the evt files both -before_SOBI and After_SOBI
% 2- in the before_SOBI file , delete the last column ( Comment Column ) -
% you can do it by copy and paste all the file to Excel 
% 3- Get rid of any trigger No (3rd column in the original file), in the before_SOBI file, which is not in
% the range of [11-17 and 101] for MSI and [31-34  and 444,555] for ANTI.
% 4- Check the number of the stim trigger matches with the number of
% response trig. Otherwise the code gives you an error 
% 5- Check for any blank lines in the evt files (last lines) and delete them 
% Run the code


before_sobi_evt='C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\6030_anti_stim_-500+1100_forSOBI_edited_2.evt';
after_sobi_evt='C:\Iman Work\CMB Projects\Yukari Project\6030_anti_stim\recon_6030_anti_stim_-500+1100_forSOBI_ACCA.evt';

% reading the evt files 

[tmu1,code1,triger1]=textread(before_sobi_evt,'%d %d %d');
[tmu2,code2,triger2,desc2]=textread(after_sobi_evt,'%d %d %d %s');
mat_new=[];
[r1,c1]=size(tmu1);
[r2,c2]=size(tmu2);

% Checking for compatibility of the evt file
if r1~=r2*2 
    disp('sth wrong in the forSOBI file');
    quit
end  


% computing the new evt file
counter=0;
for i=1:r2
    diff=tmu1(2*i)-tmu1(2*i-1);
    mat_new=[mat_new;tmu2(i),code2(i),triger2(i),desc2(i);tmu2(i)+diff, code1(2*i),triger1(2*i),desc2(i)];
end
[r3,c3]=size(mat_new);


% Writing the new evt file
fid = fopen(strcat(strrep(after_sobi_evt,'.evt',''),'_','AutoTrigerFixed.evt'),'w');
if fid == -1
   fprintf(1,'Error creating New Event File  \n');
end
fprintf(fid,'Tmu\tCode\tTriNo\tComnt\n');
for i=1:r3
            fprintf(fid, '%d\t%d\t%d\t%s\n',cell2mat(mat_new(i,1)), cell2mat(mat_new(i,2)),cell2mat(mat_new(i,3)),cell2mat(mat_new(i,4))) ;
end 

    fclose(fid);    
