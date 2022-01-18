function string = checkPhase(string)
%Checks the phase of the supplied string.
%The first bit of the string is used for the calculation.
%
%word = checkPhase(string)
%
%   Inputs:
%       string        - an array with 85 bit of data from the navigation
%                       message (a character array, must contain only '0' or
%                       '1'). 
%
%   Outputs:
%       string        - string with corrected polarity of the data bits
%                       (character array). 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Written by Darius Plausinaitis and Dennis M. Akos
%
% GLONASS modification by Jakob Almqvist
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
% $Id: checkPhase.m,v 1.1.2.4 2006/08/14 11:38:22 dpl Exp $

if string(1) == '1'
    % Data bits must be inverted
    string = invert(string);
end
