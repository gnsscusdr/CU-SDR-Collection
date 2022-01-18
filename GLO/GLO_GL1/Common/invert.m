function result = invert(data)
% Inverts the binary input-string so that 0 becomes 1 and 1 becomes 0.
%
%result = invert(data)

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Written by Darius Plausinaitis, Kristin Larson and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

% CVS record:
% $Id: invert.m,v 1.1.2.4 2006/08/14 11:38:22 dpl Exp $

% !!!!! Important comment !!!!!
%This function differs from the other versions of the invert function
%across the other SDRs. Old function was usin dec2bin matlab function to
%invert the data. Problem, in Glonass we need to invert data string of
%length 85, which is too much for dec2bin (max : 55). The use of sscanf
%takes care of the problem. If any question, please have a look at :
%https://www.mathworks.com/help/matlab/ref/bin2dec.html

dataArray = sscanf(data,'%1x')';
result = sprintf('%d',1-dataArray);
