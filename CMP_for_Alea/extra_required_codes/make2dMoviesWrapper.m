% Aim: to make a wrapper around make2dmovies() so that multiple movies can
% be launched at a time.
%
% Input: xlsfilename, the name of the xlsfile (full path). This file has
% 5 columns: (1) full name (including path) of csdmul file, 
%            (2) full name (including path) of avgmul file, 
%            (3) full name (including path) of sfp file.
%            (4) min scale value (leave them blank for [])
%            (5) max scale value (leave them blank for [])
%
% Output: the generated movies will be stored in the folder where mul files
% for the individual subject is stored.
%
%
% Example: make2dmoviesWrapper('/Users/Shared/autism/make2dmoviesWrapper_xls.xls')
%
% Written by Manish Saggar (mishu@cs.utexas.edu) on 4/15/2011

function y = make2dMoviesWrapper(xlsfilename)

    [n,t,r] = xlsread(xlsfilename);
    
    % first line in the xls file is for the header.
    for i = 2:1:size(r,1)
        [pathstr,name,ext] = fileparts(r{i,1});
        cd(pathstr);
        if isnan(r{i,5}) || isnan(r{i,6})
            make2dMovies(r{i,1}, r{i,2}, r{i,3}, r{i,4}, [], []);
        else
            make2dMovies(r{i,1}, r{i,2}, r{i,3}, r{i,4}, r{i,5}, r{i,6});
        end
        cd /Users/Shared/autism/;
    end

end