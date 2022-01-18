function L5Qcode = generateL5Qcode(PRN,settings)
% generateL2Ccode.m generates one of the GPS satellite L2C codes.
%
% L2Ccode = generateL2Ccode(PRN,flag,settings)
%
%   Inputs:
%       PRN         - PRN number of the sequence. 
%       settings    - receiver settings
%
%   Outputs:
%       L5Icode      - a vector containing the desired L5I code sequence 
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
%$Id: generateL2Ccode.m,v 1.1.2.5 2017/01/16 23:00:00 dpl Exp $


%--- Code length --------------------------------------------------------------------
CodeLength = settings.codeLength;

%--- Generate XA codes --------------------------------------------------------------

% Feedback position of exclusive OR operation for XA codes
xa_FeedbackPos = [9 10 12 13];

% Initial state of XA register: 1 for 0, and -1 for 1, to perform exclusive 
% OR operation by multiply
xa_reg = ones(1,13) * -1;

% XA codes output
XA = zeros(1,CodeLength);

% XA register will be reset to all 1 (here -1 represents 1)
reset_state = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1];

for ind = 1:CodeLength
    XA(ind) = xa_reg(end);
    if isequal(xa_reg,reset_state)
        xa_reg = ones(1,13) * -1;
    else
        % Exclusive OR operation for feedback
        feedback = prod(xa_reg(xa_FeedbackPos));
        % shift the register to right by one element
        xa_reg = circshift(xa_reg',1)';
        xa_reg(1) = feedback;        
    end
end

%--- Generate XBI codes -------------------------------------------------------------
% Initial-state table from pages 5--7 and pages 29--33 of IS-GPS-705D
l5q_init = [
    1701,   323,    5292,   2020,   5429,   7136,   1041,   5947,...
    4315,   148,    535,    1939,   5206,   5910,   3595,   5135,...
    6082,   6990,   3546,   1523,   4548,   4484,   1893,   3961,...
    7106,   5299,   4660,   276,    4389,   3783,   1591,   1601,...
    749,    1387,   1661,   3210,   708,    4226,   5604,   6375,...
    3056,   1772,   3662,   4401,   5218,   2838,   6913,   1685,...
    1194,   6963,   5001,   6694,   991,    7489,   2441,   639,...
    2097,   2498,   6470,   2399,   242,    3768,   1186,   5246,...
    4259,   5907,   3870,   3262,   7387,   3069,   2999,   7993,...
    7849,   4157,   5031,   5986,   4833,   5739,   7846,   898,...
    2022,   7446,   6404,   155,    7862,   7795,   6121,   4840,...
    6585,   429,    6020,   200,    1664,   1499,   7298,   1305,...
    7323,   7544,   4438,   2485,   3387,   7319,   1853,   5781,...
    1874,   7555,   2132,   6441,   6722,   1192,   2588,   2188,...
    297,    1540,   4138,   5231,   4789,   659,    871,    6837,...
    1393,   7383,   611,    4920,   5416,   1611,   2474,   118,...
    1382,   1092,   7950,   7223,   1769,   4721,   1252,   5147,...
    2165,   7897,   4054,   3498,   6571,   2858,   8126,   7017,...
    1901,   181,    1114,   5195,   7479,   4186,   3904,   7128,...
    1396,   4513,   5967,   2580,   2575,   7961,   2598,   4508,...
    2090,   3685,   7748,   684,    913,    5558,   2894,   5858,...
    6432,   3813,   3573,   7523,   5280,   3376,   7424,   2918,...
    5793,   1747,   7079,   2921,   2490,   4119,   3373,   977,...
    681,    4273,   5419,   5626,   1266,   5804,   2414,   6444,...
    4757,   427,    5452,   5182,   6606,   6531,   4268,   3115,...
    6835,   862,    4856,   2765,   37,     1943,   7977,   2512,...
    4451,   4071 ];

% Feedback position of exclusive OR operation for XBI codes
xbq_FeedbackPos = [1 3 4 6 7 8 12 13];

% Initial state of XBI register: 1 for 0, and -1 for 1, to perform exclusive 
% OR operation by multiply
xbq_reg = ones(1,13) * -1;

% XA codes output
XBQ = zeros(1,CodeLength);

% XBI code Advance(Chips) for the start of the XBI code
resetPos = l5q_init(PRN);

% XBI initial states
for ind = 1:resetPos
    % Exclusive OR operation for feedback
    feedback = prod(xbq_reg(xbq_FeedbackPos));
    % shift the register to right by one element
    xbq_reg = circshift(xbq_reg',1)';
    xbq_reg(1) = feedback; 
end

% XBQ output
for ind = 1:CodeLength
    XBQ(ind) = xbq_reg(end);
    % Exclusive OR operation for feedback
    feedback = prod(xbq_reg(xbq_FeedbackPos));
    % shift the register to right by one element
    xbq_reg = circshift(xbq_reg',1)';
    xbq_reg(1) = feedback; 
end

%---  L5I PRN code output -----------------------------------------------------------
L5Qcode = XBQ .* XA;
