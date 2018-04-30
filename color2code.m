function [code] =  color2code(color)

if  (strcmp(color,'red') == 1)
    code = 2;
elseif (strcmp(color,'blue') == 1)
    code = 1;
elseif (strcmp(color,'green') == 1)
    code = 16;
elseif (strcmp(color,'yellow') == 1)
    code = 8;
elseif (strcmp(color,'none') == 1)
    code = 500;
else    
    disp('Response code not matched')
end

end