function Secondary = generatePilot2ndCodes(PRN)
% This function generates Galileo E1B code in bipolar format (-1, +1)
% PRNs 1 to 50 (no care is taken if PRN number is > 50)
% The memory codes are stored in the file gale1bcode.dat
%
% E1Bcode = generateE1Bcode(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%
%   Outputs:
%       E5BI        - a vector containing the desired E5BI code sequence
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



wp_pilot = [269 	1889; 1448 	1268; 1028 	1593; 1324 	1186; ... 
822 	1239; 5 	1930; 155 	176; 458 	1696; 310 	26; ... 
959 	1344; 1238 	1271; 1180 	1182; 1288 	1381; 334 	1604; ... 
885 	1333; 1362 	1185; 181 	31; 1648 	704; 838 	1190; ... 
313 	1646; 750 	1385; 225 	113; 1477 	860; 309 	1656; ... 
108 	1921; 1457 	1173; 149 	1928; 322 	57; 271 	150; ... 
576 	1214; 1103 	1148; 450 	1458; 399 	1519; 241 	1635; ... 
1045 	1257; 164 	1687; 513 	1382; 687 	1514; 422 	1; ... 
303 	1583; 324 	1806; 495 	1664; 725 	1338; 780 	1111; ... 
367 	1706; 882 	1543; 631 	1813; 37 	228; 647 	2871; ... 
1043 	2884; 24 	1823; 120 	75; 134 	11; 136 	63; ... 
158 	1937; 214 	22; 335 	1768; 340 	1526; 661 	1402; ... 
889 	1445; 929 	1680; 1002 	1290; 1149 	1245];
        
        

N = 3607;

legendre = zeros(1,N);
for ind = 1:N-1
    legendre(ind+1) = JacobiSymbol(ind,N);
end

legendre(legendre==-1) = 0;

% The length of the secondary code of B1C is 1800
Secondary = zeros(1,1800);

p = wp_pilot(PRN,2);
w = wp_pilot(PRN,1);

for ind=0:1799

    k = mod((ind+p-1), N);
    
    Secondary(ind+1) = xor( legendre(k+1), legendre(mod(k+w,N)+1) );
end


% Change to bipolar format (-1, +1)
Secondary = 1-2*Secondary;


%% Test the generation
% % Secondary = Secondary(end-23:end);
% Secondary = (1-Secondary)/2;
% dicmal = zeros(1,8);
% for ind = 1:8
%     a = Secondary((ind-1)*3+1:(ind-1)*3+3);
%     dicmal(ind)  = a(1)*4 + a(2)*2 + a(3);
% end        
%         
%  dicmal






