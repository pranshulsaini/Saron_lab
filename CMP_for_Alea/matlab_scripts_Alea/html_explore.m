
%add path for tif files
dir = 'C:\Users\plsaini\Box Sync\Stroop\Temp\matlab_scripts\TIF_files\';

% print the HTML file
fid = fopen('ArtifactsHTML.html','w');
fprintf(fid, '\n<HTML>\n<TITLE>%s</TITLE>\n<BODY>\n <form action="http://www.rebol.com/cgi-bin/test-cgi.cgi" method="POST">\n<TABLE border=1>\n', strcat('try','_','STR'));
fprintf(fid, '<input type="text" name="filename", value=%s>\n','try_STR');
fprintf(fid, '<input type="text" name="nsources", value=%s>\n',num2str(20));

%store an Excel file that lists the sources in their order of appearance. To be uesd for Voting%%
sourceorder = [];
%fileid = fopen(strcat(strrep(HTMLtitle,'.dat',''),'_voteforrecon.xls'),'w');

if fid == -1
    fprintf(1,'Error creating Voting file in Excel\n');
end


emg2 = 4;    
for i = 1:emg2
    if emg2 > 3
        fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm"><input type="radio" name="%d" value="emg" ><input type="radio" name="%d" value="Rnorm">', num2str(i),2*i,3*i,4*i);
    else
        fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm" ><input type="radio" name="%d" value="emg"><input type="radio" name="%d" value="Rnorm">', num2str(10*i),20*i,30*i,40*i);
    end
    fprintf(fid, '\n<IMG SRC=''%s''></TD></TR>', strcat(dir,'emg_',num2str(i),'.tif')); 
    sourceorder = [sourceorder i];
end

eog2 = 5;
for i = 1:eog2
    if eog2 > 2
        fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2><input type="radio" name="%d" value="norm"><input type="radio" name="%d" value="eog" ><input type="radio" name="%d" value="Rnorm">', num2str(i),2*i,3*i,4*i);        
    else
        fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2><input type="radio" name="%d" value="norm" ><input type="radio" name="%d" value="eog"><input type="radio" name="%d" value="Rnorm">', num2str(10*i),20*i,30*i,40*i);        
    end
    fprintf(fid, '\n<IMG SRC=''%s''></TD></TR>', strcat(dir,'eog_',num2str(i),'.tif'));  
    sourceorder = [sourceorder i];
end

normal2 = 6;
for i = 1:normal2
    if normal2 > 3
        fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value=norm><input type="radio" name="%d" value="eog" ><input type="radio" name="%d" value="Rnorm">', num2str(i),2*i,3*i,4*i);        
    else
        fprintf(fid, '\n<TR><TD><P><H2 align=center> Source = %s </H2> <input type="radio" name="%d" value="norm" ><input type="radio" name="%d" value="emg"><input type="radio" name="%d" value="Rnorm">', num2str(10*i),20*i,30*i,40*i);       
    end
    fprintf(fid, '\n<IMG SRC=''%s''></TD></TR>', strcat(dir,'normal_',num2str(i),'.tif'));  
    sourceorder = [sourceorder i];
end
fprintf(fid, '\n</TABLE><input type="submit" value="Vote"></BODY></HTML>');

fclose(fid);
xlswrite(strcat('try','_voteforrecon.xlsx'),sourceorder');