function update_sfp(sfpfile)

fid = fopen(sfpfile);
s = textscan(fid, '%s %f %f %f');
fclose(fid);

y = 0;
for i = 1: size(s{1,1},1)
    x = s{1,1}(i);
    if length(x{1,1}) > 5
        y = strcmp(x{1,1}(4:6), 'AVR');
    end
    if y
        break;
    end
end

n_headerlines = i - 1;

fid = fopen(sfpfile);
MyText = textscan(fid,'%s %f %f %f','headerlines',n_headerlines);
fclose(fid);

n_rows = size(MyText{1,1},1);

col1 = strings(n_rows, 1 );
col1(:) = MyText{1,1}(:);
col2 = MyText{1,2}(:);
col3 = MyText{1,3}(:);
col4 = MyText{1,4}(:);


[tmp1, tmp2, tmp3] = fileparts(sfpfile); 
filename = strcat(tmp1,'\',tmp2,'_updated',tmp3);
fid = fopen(filename,'w');

for i = 1:n_rows
    fprintf(fid, '%s\t%f\t%f\t%f\n', col1(i), col2(i), col3(i), col4(i));
end
fclose(fid);

end
