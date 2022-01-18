function E5aI = generateE5aIcode(PRN,flag)
% This function generates Galileo E5aI primary/tiered code in bipolar
% format (-1, +1) for PRNs 1 to 50 (no care is taken if PRN number is > 50)
%
% E5aI = generateE5aIcode(PRN,flag)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%       flag        - 1: primary code; 2: tiered code.
%
%   Outputs:
%       E5aI        - a vector containing the E5aI primary/tiered code
%                   sequence (chips).

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

% Start values for Register 1 (see Galileo-OS-SIS-ICD) are all 1s.
Register1 = ones(1,14);

% Start values for Register 2 (see Galileo-OS-SIS-ICD)
e5ai_init = [ ...
    '30305',    '14234',    '27213',    '20577', ...
    '23312',    '33463',    '15614',    '12537', ...
    '01527',    '30236',    '27344',    '07272', ...
    '36377',    '17046',    '06434',    '15405', ...
    '24252',    '11631',    '24776',    '00630', ...
    '11560',    '17272',    '27445',    '31702', ...
    '13012',    '14401',    '34727',    '22627', ...
    '30623',    '27256',    '01520',    '14211', ...
    '31465',    '22164',    '33516',    '02737', ...
    '21316',    '35425',    '35633',    '24655', ...
    '14054',    '27027',    '06604',    '31455', ...
    '34465',    '25273',    '20763',    '31721', ...
    '17312',    '13277'];

% Feedback Tap coefficients for Register 1 & Register 2
Feedback_Reg1 = '40503';
Feedback_Reg2 = '50661';

% Covert to decimal format. The ASCII code for '0' is 48,
% and for '1' is 49
taps1_coef = dec2bin(base2dec(Feedback_Reg1,8))-48;
taps2_coef = dec2bin(base2dec(Feedback_Reg2,8))-48;

% Take the 14 significant bits for feedback Tap coefficients from higher
% bits to lower bits: 14th, 13th, 12th ... 1st
taps1_coef = taps1_coef(1:14);
taps2_coef = taps2_coef(1:14);

% Start values for Register 2
StartValues = e5ai_init(5*(PRN-1)+1:5*PRN);
StartValues = dec2bin(base2dec(StartValues,8))-48;

% Consider when the first bit in the start values for register 2 is 0:
% dec2bin() and base2dec() omit the higher bit of 0
Register2 = zeros(1,14);
Register2(end-length(StartValues)+1:end) = StartValues;

% E5aI primary code
Pri_E5aI = zeros(1,10230);
for ind = 1:length(Pri_E5aI)
    
    % bipolar format output of the code generator
    RegOut1 = Register1 .* taps1_coef;
    RegOut2 = Register2 .* taps2_coef;
    Pri_E5aI(ind) = (1 - 2* RegOut1(1)) * (1 - 2* RegOut2(1));
    
    % Exclusive OR operation for feedback: change to bipolar format (-1,
    % +1) to perform exclusive OR operation by multiply: (0-->1; 1-->-1)
    feedback1 = prod(1- 2* RegOut1);
    feedback2 = prod(1- 2* RegOut2);
    
    % Re-convert to original format (1-->0; -1-->1)
    if feedback1 == 1
        feedback1 = 0;
    elseif feedback1 == -1
        feedback1 = 1;
    end
    if feedback2 == 1
        feedback2 = 0;
    elseif feedback2 == -1
        feedback2 = 1;
    end
    
    % shift the register to left ( to higher bit) by one element
    Register1 = circshift(Register1',-1)';
    Register2 = circshift(Register2',-1)';
    % Shift the feedbacked value into the lowerest bit
    Register1(end) = feedback1;
    Register2(end) = feedback2;
end

% ----- E5aI primary or tiered code generation ----------------------------
if flag ==1
    E5aI = Pri_E5aI;
elseif  flag ==2
    % Secondary code is "842E9"
    SecondaryCode = [-1 1 1 1   1 -1 1 1   1 1 -1 1   -1 -1 -1 1   -1 1 1 -1];
    
    % -- Combine the promary and secondary codes to form the tiered code --
    E5aI = zeros(1,10230*20);
    for ind = 1:20
        % 20ms E5aI code
        E5aI( (ind-1)*10230+1 : ind*10230 ) = Pri_E5aI * SecondaryCode(ind);
    end
end
