function CMcode = generateCMcode(PRN,settings)
% generateL2Ccode.m generates one of the GPS satellite L2C return-zero codes .
%
% L2Ccode = generateL2Ccode(PRN,flag,settings)
%
%   Inputs:
%       PRN         - PRN number of the sequence. 
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
%--------------------------------------------------------------------------

%CVS record:
%$Id: generateL2Ccode.m,v 1.1.2.5 2016/08/23 22:00:00 dpl Exp $

%--- Initialize related parameters ----------------------------------------
RegPos = [4 7 9 12 15 17 19 22 23 24 25];

%--- Initial-state table from pages 9--11 and pages 62--63 of IS-GPS-200H--
% PRNs 64-158 do not exist for L2 CM-/L2 CL-code.
l2cm_init = [...
   ... PRNs 1-37 for GPS satellites
   742417664,   756014035,   002747144,   066265724,...
   601403471,   703232733,   124510070,   617316361,...
   047541621,   733031046,   713512145,   024437606,...
   021264003,   230655351,   001314400,   222021506,...
   540264026,   205521705,   064022144,   120161274,...
   044023533,   724744327,   045743577,   741201660,...
   700274134,   010247261,   713433445,   737324162,...
   311627434,   710452007,   722462133,   050172213,...
   500653703,   755077436,   136717361,   756675453,...
   435506112,...
   ...PRNs 38-63 are required per this Table if a manufacturer chooses to 
   ...include these PRNs in their receiver design
   771353753,   226107701,   022025110,   402466344,...
   752566114,   702011164,   041216771,   047457275,...
   266333164,   713167356,   060546335,   355173035,...
   617201036,   157465571,   767360553,   023127030,...
   431343777,   747317317,   045706125,   002744276,...
   060036467,   217744147,   603340174,   326616775,...
   063240065,   111460621,...
   ... PRN 159-210 are reserved for other GNSS applications.
   604055104,   157065232,   013305707,   603552017,...
   230461355,   603653437,   652346475,   743107103,...
   401521277,   167335110,   014013575,   362051132,...
   617753265,   216363634,   755561123,   365304033,...
   625025543,   054420334,   415473671,   662364360,...
   373446602,   417564100,   000526452,   226631300,...
   113752074,   706134401,   041352546,   664630154,...
   276524255,   714720530,   714051771,   044526647,...
   207164322,   262120161,   204244652,   202133131,...
   714351204,   657127260,   130567507,   670517677,...
   607275514,   045413633,   212645405,   613700455,...
   706202440,   705056276,   020373522,   746013617,...
   132720621,   434015513,   566721727,   140633660];

%--- Select initial states for PRN ----------------------------------------
if (PRN >= 1 && PRN <= 63)
    shiftPos = PRN;
elseif(PRN >= 159 && PRN <= 210)
    shiftPos = PRN - 95;
elseif ( PRN < 1 || PRN > 210 || (PRN >= 64 && PRN <= 158) )
    disp('  PRN does not exist! ');
    return;
end

%--- Code length and initial state ----------------------------------------
CodeLength = settings.codeLength;
code_init = l2cm_init(shiftPos);
   
%--- Covert octal?format to binary bit and fill in shift register ---------
reg = dec2bin(oct2dec(code_init)) - 48;
reg = [zeros(1,27-length(reg)) reg];

%--- Convert 0 to 1 and  1 to -1 ------------------------------------------
reg(find(reg)) = -1;  %#ok<*FNDSB>
reg(find(reg == 0)) = 1;

%--- CM code generating ---------------------------------------------------
CMcode = zeros(1,CodeLength);
for index = 1:CodeLength
    CMcode(index) = reg(end);
    reg = circshift(reg',1)';
    reg(RegPos) = reg(RegPos) * CMcode(index); 
end

% --- return-zero CM codes ------------------------------------------------
CMcode = [CMcode;zeros(1,length(CMcode))];
CMcode = reshape(CMcode,1,[]);

