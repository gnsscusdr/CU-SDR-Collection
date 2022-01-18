function E5aQ = generateE5aQcode(PRN,flag)
% This function generates Galileo E5aQ primary/tiered code in bipolar 
% format (-1, +1) for PRNs 1 to 50 (no care is taken if PRN number is > 50)
%
% E5aQ = generateE5aQcode(PRN)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%       flag        - 1: primary code; 2: tiered code.
%
%   Outputs:
%       E5BQ        - a vector containing the E5aQ primary code sequence
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
%$Id: generateE5BIcode.m,v 1.1.2.5 2017/04/24 22:00:00 dpl Exp $

%% Prepare for related parameters ===================================
% Start values for Register 1 (see Galileo-OS-SIS-ICD) are all 1s.
Register1 = ones(1,14);

% Initial-state table from Galileo-OS-SIS-ICD (2016) for Register 2
e5aq_init = [ ...
   '25652',    '05142',    '24723',    '31751', ...
   '27366',    '24660',    '33655',    '27450', ...
   '07626',    '01705',    '12717',    '32122', ...
   '16075',    '16644',    '37556',    '02477', ...
   '02265',    '06430',    '25046',    '12735', ...
   '04262',    '11230',    '00037',    '06137', ...
   '04312',    '20606',    '11162',    '22252', ...
   '30533',    '24614',    '07767',    '32705', ...
   '05052',    '27553',    '03711',    '02041', ...
   '34775',    '05274',    '37356',    '16205', ...
   '36270',    '06600',    '26773',    '17375', ...
   '35267',    '36255',    '12044',    '26442', ...
   '21621',    '25411'];

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
StartValues = e5aq_init(5*(PRN-1)+1:5*PRN);
StartValues = dec2bin(base2dec(StartValues,8))-48;

% Consider when the first bit in the start values for register 2 is 0:
% dec2bin() and base2dec() omit the higher bit of 0
Register2 = zeros(1,14);
Register2(end-length(StartValues)+1:end) = StartValues;

%% Generating the E5aQ code =========================================
Pri_E5aQ = zeros(1,10230);
for ind = 1:length(Pri_E5aQ)
    % Register output mutiplied by the coefficeients
    RegOut1 = Register1 .* taps1_coef;
    RegOut2 = Register2 .* taps2_coef;
    
    % E5BQ code sequence output (bipolar format)
    Pri_E5aQ(ind) = (1 - 2* RegOut1(1)) * (1 - 2* RegOut2(1));
    
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

% ----- E5aQ primary or tiered code generation ----------------------------
if flag ==1
    E5aQ = Pri_E5aQ;
elseif  flag ==2
    % Secondary code
    SecondaryCode = generateE5aQ_secondary(PRN);
    
    % -- Combine the promary and secondary codes to form the tiered code --
    E5aQ = zeros(1,10230*100);
    for ind = 1:100
        % 100ms E5aQ code
        E5aQ( (ind-1)*10230+1 : ind*10230 ) = Pri_E5aQ * SecondaryCode(ind);
    end
end