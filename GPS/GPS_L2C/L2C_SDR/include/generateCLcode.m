function CLcode = generateCLcode(PRN,settings)
% generateL2Ccode.m generates one of the GPS satellite L2C codes.
%
% L2Ccode = generateL2Ccode(PRN,flag,settings)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%       flag        - 'CM' or 'cm' for L2c CMcode; 'CL' or 'cl' for L2c CLcode; 
%       settings    - receiver settings
%
%   Outputs:
%       L2Ccode      - a vector containing the desired L2C code sequence 
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
%------------------------------------------------------------------------------------

%CVS record:
%$Id: generateL2Ccode.m,v 1.1.2.5 2016/08/23 22:00:00 dpl Exp $

%--- Initialize related parameters --------------------------------------------------
RegPos = [4 7 9 12 15 17 19 22 23 24 25];

%--- Initial-state table from pages 9--11 and pages 62--63 of IS-GPS-200H ---
% PRNs 64-158 do not exist for L2 CM-/L2 CL-code.
l2cl_init = [...
   ... PRNs 1-37 for GPS satellites
   624145772,   506610362,   220360016,   710406104,...
   001143345,   053023326,   652521276,   206124777,...
   015563374,   561522076,   023163525,   117776450,...
   606516355,   003037343,   046515565,   671511621,...
   605402220,   002576207,   525163451,   266527765,...
   006760703,   501474556,   743747443,   615534726,...
   763621420,   720727474,   700521043,   222567263,...
   132765304,   746332245,   102300466,   255231716,...
   437661701,   717047302,   222614207,   561123307,...
   240713073,...
   ...For PRNs 38-63
   101232630,   132525726,   315216367,   377046065,...
   655351360,   435776513,   744242321,   024346717,...
   562646415,   731455342,   723352536,   000013134,...
   011566642,   475432222,   463506741,   617127534,...
   026050332,   733774235,   751477772,   417631550,...
   052247456,   560404163,   417751005,   004302173,...
   715005045,   001154457,...
   ... For PRN 159-210
   605253024,   063314262,   066073422,   737276117,...
   737243704,   067557532,   227354537,   704765502,...
   044746712,   720535263,   733541364,   270060042,...
   737176640,   133776704,   005645427,   704321074,...
   137740372,   056375464,   704374004,   216320123,...
   011322115,   761050112,   725304036,   721320336,...
   443462103,   510466244,   745522652,   373417061,...
   225526762,   047614504,   034730440,   453073141,...
   533654510,   377016461,   235525312,   507056307,...
   221720061,   520470122,   603764120,   145604016,...
   051237167,   033326347,   534627074,   645230164,...
   000171400,   022715417,   135471311,   137422057,...
   714426456,   640724672,   501254540,   513322453];

if (PRN >= 1 && PRN <= 63)
    shiftPos = PRN;
elseif(PRN >= 159 && PRN <= 210)
    shiftPos = PRN - 95;
elseif ( PRN < 1 || PRN > 210 || (PRN >= 64 && PRN <= 158) )
    disp('  PRN does not exist! ');
    return;
end

% 
CodeLength = settings.CLCodeLength;
code_init = l2cl_init(shiftPos);
      
reg = dec2bin(oct2dec(code_init)) - 48;
reg = [zeros(1,27-length(reg)) reg];

%--- Convert 0 to 1 and  1 to -1 ---
reg(find(reg)) = -1;  %#ok<*FNDSB>
reg(find(reg == 0)) = 1;

CLcode = zeros(1,CodeLength);
for index = 1:CodeLength
    CLcode(index) = reg(end);
    reg = circshift(reg',1)';
    reg(RegPos) = reg(RegPos) * CLcode(index); 
end

% --- return-zero CL codes ------------------------------------------------
CLcode = [zeros(1,length(CLcode));CLcode];
CLcode = reshape(CLcode,1,[]);



