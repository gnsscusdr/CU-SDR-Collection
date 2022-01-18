function CAcode = generateCAcode53(PRN)
% generateCAcode.m generates one of the 37 Beidou satellite B1I codes.
%
% CAcode = generateCAcodeBDS(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       CAcode      - a vector containing the desired C/A code sequence 
%                   (chips).  

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Daehee Won, Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

%CVS record:
%$Id: generateCAcode.m,v 1.1.2.5 2006/08/14 11:38:22 dpl Exp $

%--------------------------------------------------------------------------
% Modified for Beidou by Daehee Won
% Final update: Sep. 26, 2013

%--- Generate G1 code
% Initialize G1 output
g1 = zeros(1, 2046);    % Chip length of Beidou is 2046
% Load shift register (11 bits)
reg = -1*[-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1];

% Generate G1 signal chips based on the G1 generator polynomial
for i=1:2046
    g1(i)     = reg(11);
    saveBit   = reg(1)*reg(7)*reg(8)*reg(9)*reg(10)*reg(11);
    reg(2:11) = reg(1:10);
    reg(1)    = saveBit;
end

%--- Generate G2 code
% Initialize G2 output
g2 = zeros(1, 2046);    % Chip length of Beidou is 2046
% Load shift register (11 bits)
reg = -1*[-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1];
% Phase assignment of G2 sequence

g2s1 = [1, 1, 1, 1, 1, 1, 1, 1, 2, 3, ...
        3, 3, 3, 3, 3, 3, 4, 4, 4, 4, ...
        4, 4, 5, 5, 5, 5, 5, 6, 6, 6, ...
        6, 8, 8, 8, 9, 9, 10,... 
        repmat([1], 1 , 16),...
        repmat([2], 1 , 3), ...
        repmat([3], 1 , 2)];
    
g2s2 = [ 3, 4, 5, 6, 8, 9,10,11, 7, 4, ...
         5, 6, 8, 9,10,11, 5, 6, 8, 9, ...
        10,11, 6, 8, 9,10,11, 8, 9,10, ...
        11, 9,10,11,10,11,11, ...
        2, ...
        repmat([3], 1, 5),...
        repmat([4], 1, 2),...
        repmat([5], 1, 4),...
        6, 8, 9, 9, 3, 5, 7, 4, 4];
        
g2s3 = [7 4 6 8 10 11 5 9 6 ...
        8 10 11 9 9 10 11 7 ...
        7 9 5 9];

% Add new shift register for PRN > 37
if PRN>37
    
    for i=1:2046
        g2(i)       = reg(g2s1(PRN))*reg(g2s2(PRN))*reg(g2s3(PRN-37));
        saveBit     = reg(1)*reg(2)*reg(3)*reg(4)*reg(5)*reg(8)*reg(9)*reg(11);
        reg(2:11)   = reg(1:10);
        reg(1)      = saveBit;
    end
    
else
% Generate G2 signal chips based on the G2 generator polynomial
    for i=1:2046
        g2(i)       = reg(g2s1(PRN))*reg(g2s2(PRN));
        saveBit     = reg(1)*reg(2)*reg(3)*reg(4)*reg(5)*reg(8)*reg(9)*reg(11);
        reg(2:11)   = reg(1:10);
        reg(1)      = saveBit;
    end
end
%--- Form single sample C/A code by multiplying G1 and G2 -----------------
CAcode = -(g1 .* g2);
