function E1Bcode = generateE1Bcode(PRN)
% This function generates Galileo E1B BOC(1,1) code in bipolar 
% format (-1, +1) PRNs 1 to 50 (no care is taken if PRN number is > 50)
% The memory codes are stored in the file gale1bcode.dat
%
% E1Bcode = generateE1Bcode(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       E1Bcode      - a vector containing the desired E1B code sequence
%                   (chips).

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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
%$Id: generateE1Bcode.m,v 1.1.2.5 2017/04/24 22:00:00 dpl Exp $

%--- Generate the primary code --------------------------------------------
% primary code length
PRIMARY_CODE_LENGTH = 4092;

% Number of Galileo E1B PRN codes
NUM_GAL_PRNS = 50;

% File for the memery code
fid = fopen('E1b.dat');

% Read the code sequenrences out
allprncode = fscanf(fid,'%d', PRIMARY_CODE_LENGTH*NUM_GAL_PRNS);

% Take the code sequenrences for the specified PRN number
e1b_primary_code= allprncode(...
                   (PRN-1)*PRIMARY_CODE_LENGTH+1:PRN*PRIMARY_CODE_LENGTH);

% Change to bipolar format (-1, +1)
e1bCodeRaw = 1-2*e1b_primary_code;

%--- Add BOC(1,1) subcarrier  ---------------------------------------------
E1Bcode = zeros(1, length(e1bCodeRaw)*2);
jj = 1;
for ii = 1:2:length(E1Bcode)-1
    E1Bcode(ii)   =  e1bCodeRaw(jj);
    E1Bcode(ii+1) = - e1bCodeRaw(jj);
    jj = jj + 1;
end
% Close the fid
fclose(fid);




