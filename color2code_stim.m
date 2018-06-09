function [code] =  color2code(color)
% this is for responses


if  (strcmp(color,'red') == 1)
    code = 6;
elseif (strcmp(color,'blue') == 1)
    code = 5;
elseif (strcmp(color,'green') == 1)
    code = 20;
elseif (strcmp(color,'yellow') == 1)
    code = 12;
else    
    disp('Response code not matched')
end

end