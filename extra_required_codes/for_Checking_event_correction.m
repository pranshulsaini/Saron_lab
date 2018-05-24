for i=2:r2
    
    a1=cell2mat(triger1(i*2-2)); a1=str2num(a1);
    a2=cell2mat(triger2(i)); a2=str2num(a2);
    a=a1-a2;
    b=num2str(a); c=num2str(i);
    d=strcat(c,'>>', b);
    disp(d)
    
    
end
