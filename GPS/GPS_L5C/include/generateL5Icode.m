function L5Icode = generateL5Icode(PRN,settings)
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
l5i_init = [ ...
    266,    365,    804,    1138,   1509,   1559,   1756,   2084,  ...
    2170,   2303,   2527,   2687,   2930,   3471,   3940,   4132,  ...
    4332,   4924,   5343,   5443,   5641,   5816,   5898,   5918, ...
    5955,   6243,   6345,   6477,   6518,   6875,   7168,   7187,...
    7329,   7577,   7720,   7777,   8057,   5358,   3550,   3412,...
    819,    4608,   3698,   962,    3001,   4441,   4937,   3717, ...
    4730,   7291,   2279,   7613,   5723,   7030,   1475,   2593,...
    2904,   2056,   2757,   3756,   6205,   5053,   6437,   7789,...
    2311,   7432,   5155,   1593,   5841,   5014,   1545,   3016,...
    4875,   2119,   229,    7634,   1406,   4506,   1819,   7580,...
    5446,   6053,   7958,   5267,   2956,   3544,   1277,   2996,...
    1758,   3360,   2718,   3754,   7440,   2781,   6756,   7314,...
    208,    5252,   696,    527,    1399,   5879,   6868,   217,...
    7681,   3788,   1337,   2424,   4243,   5686,   1955,   4791,...
    492,    1518,   6566,   5349,   506,    113,    1953,   2797,...
    934,    3023,   3632,   1330,   4909,   4867,   1183,   3990,...
    6217,   1224,   1733,   2319,   3928,   2380,   841,    5049,...
    7027,   1197,   7208,   8000,   152,    6762,   3745,   4723,...
    5502,   4796,   123,    8142,   5091,   7875,   330,    5272,...
    4912,   374,    2045,   6616,   6321,   7605,   2570,   2419,...
    1234,   1922,   4317,   110,    825,    958,    1089,   7813,...
    6058,   7703,   6702,   1714,   6371,   2281,   1986,   6282,...
    3201,   3760,   1056,   6233,   1150,   2823,   6250,   645,...
    2401,   1639,   2946,   7091,   923,    7045,   6493,   1706,...
    5836,   926,    6086,   950,    5905,   3240,   6675,   3197,...
    1555,   3589,   4555,   5671,   6948,   4664,   2086,   5950, ...
    5521,   1515 ];

% Feedback position of exclusive OR operation for XBI codes
xbi_FeedbackPos = [1 3 4 6 7 8 12 13];

% Initial state of XBI register: 1 for 0, and -1 for 1, to perform exclusive 
% OR operation by multiply
xbi_reg = ones(1,13) * -1;

% XA codes output
XBI = zeros(1,CodeLength);

% XBI code Advance(Chips) for the start of the XBI code
resetPos = l5i_init(PRN);

% XBI initial states
for ind = 1:resetPos
    % Exclusive OR operation for feedback
    feedback = prod(xbi_reg(xbi_FeedbackPos));
    % shift the register to right by one element
    xbi_reg = circshift(xbi_reg',1)';
    xbi_reg(1) = feedback; 
end

% XBI output
for ind = 1:CodeLength
    XBI(ind) = xbi_reg(end);
    % Exclusive OR operation for feedback
    feedback = prod(xbi_reg(xbi_FeedbackPos));
    % shift the register to right by one element
    xbi_reg = circshift(xbi_reg',1)';
    xbi_reg(1) = feedback; 
end

%---  L5I PRN code output -----------------------------------------------------------
L5Icode = XBI .* XA;
