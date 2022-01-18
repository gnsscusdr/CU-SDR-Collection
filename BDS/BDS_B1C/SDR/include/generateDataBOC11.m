function B1CData = generateDataBOC11(settings,PRN)
% This function generates BDS-3 B1C BOC(1,1) code in bipolar 
% format (-1, +1). PRNs 1 to 63 (no care is taken if PRN number is > 63)
%
% B1CData = generateDataBOC11(settings,PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       B1CData     - a vector containing the desired B1C code sequence
%                   (chips).

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos

% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and 
% implementation of an open-source BDS-3 B1C/B2a SDR receiver. 
% GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
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
%$Id: generateE5BIcode.m,v 1.1.2.5 2017/11/27 22:00:00 dpl Exp $

% Matrix of phase difference(w) and truncation point (p) for B1C 
% data primary code
wp_data = [2678   699;    4802   694;    958    7318;   859    2127;...
      3843   715;    2232   6682;   124    7850;   4352   5495;...
      1816   1162;   1126   7682;   1860   6792;   4800   9973;...
      2267   6596;   424    2092;   4192   19;     4333   10151;...
      2656   6297;   4148   5766;   243    2359;   1330   7136;...
      1593   1706;   1470   2128;   882    6827;   3202   693;...
      5095   9729;   2546   1620;   1733   6805;   4795   534;...
      4577   712;    1627   1929;   3638   5355;   2553   6139;...
      3646   6339;   1087   1470;   1843   6867;   216    7851;...
      2245   1162;   726    7659;   1966   1156;   670    2672;...
      4130   6043;   53     2862;   4830   180;    182    2663;...
      2181   6940;   2006   1645;   1080   1582;   2288   951;...
      2027   6878;   271    7701;   915    1823;   497    2391;...
      139    2606;   3693   822;    2054   6403;   4342   239;...
      3342   442;    2592   6769;   1007   2560;   310    2502;...
      4203   5072;   455    7268;   4318   341];
     
%  ---- Compute the Jacobi symbols ----------------------------------------
N = 10243;
legendre = zeros(1,N);
for ind = 1:N-1
    legendre(ind+1) = JacobiSymbol(ind,N);
end

% ---- Generate the B1C primary code --------------------------------------
legendre(legendre==-1) = 0;
Primary  = zeros(1,settings.codeLength);

% Phase difference and truncation point
p = wp_data(PRN,2);
w = wp_data(PRN,1);

% Generate the B1C primary code
for ind=0:10229
    k = mod((ind+p-1), N);
    Primary (ind+1) = xor( legendre(k+1), legendre(mod(k+w,N)+1) );
end

% Change to bipolar format (-1, +1)
Primary  = 1-2*Primary ;

% ---- Add BOC(1,1) subcarrier  -------------------------------------------
B1CData = zeros(1, length(Primary)*2);
jj = 1;
for ii = 1:2:length(B1CData)-1
    B1CData(ii)   =  -Primary(jj);
    B1CData(ii+1) = Primary(jj);
    jj = jj + 1;
end