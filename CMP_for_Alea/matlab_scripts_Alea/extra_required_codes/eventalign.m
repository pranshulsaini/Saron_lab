% eventalign - function called by pop_importevent() to find the best 
%              sampling rate ratio to align 2 arrays of data events.
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, Dec 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 9 Feb 2002, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: eventalign.m,v $
% Revision 1.2  2007/02/27 02:41:08  toby
% corrected function call, courtesy of Ronny Lindner
%
% Revision 1.1  2003/12/11 02:26:22  arno
% Initial revision
%

function mindiff = eventalign(factor, a, b)
    
       diffarray = abs(factor*a-b)';
       [allmins poss] = min(diffarray);
       mindiff = mean(allmins);
       
