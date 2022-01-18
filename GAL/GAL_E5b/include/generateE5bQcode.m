function E5bQ = generateE5bQcode(PRN,flag)
% This function generates Galileo E5bQ primary/tiered code in bipolar 
% format (-1, +1) for PRNs 1 to 50 (no care is taken if PRN number is > 50)
%
% E5BQ = generateE5bQcode(PRN,flag)
%
%   Inputs:
%       PRN         - PRN number of the sequence.
%       flag        - 1: primary code; 2: tiered code.
%
%   Outputs:
%       E5BQ        - a vector containing the E5bQ primary code sequence
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
e5bq_init2 = [...
    '03331',     '06143',     '25322',     '23371',...
    '00413',     '36235',     '17750',     '04745',...
    '13005',     '37140',     '30155',     '20237',...
    '03461',     '31662',     '27146',     '05547',...
    '02456',     '30013',     '00322',     '10761',...
    '26767',     '36004',     '30713',     '07662',...
    '21610',     '20134',     '11262',     '10706',...
    '34143',     '11051',     '25460',     '17665',...
    '32354',     '21230',     '20146',     '11362',...
    '37246',     '16344',     '15034',     '25471',...
    '25646',     '22157',     '04336',     '16356',...
    '04075',     '02626',     '11706',     '37011',...
    '27041',     '31024'];

% Feedback Tap coefficients for Register 1 & Register 2
Feedback_Reg1 = '64021';
Feedback_Reg2 = '43143';

% Covert to decimal format. The ASCII code for '0' is 48, 
% and for '1' is 49
taps1_coef = dec2bin(base2dec(Feedback_Reg1,8))-48;
taps2_coef = dec2bin(base2dec(Feedback_Reg2,8))-48;

% Take the 14 significant bits for feedback Tap coefficients from higher
% bits to lower bits: 14th, 13th, 12th ... 1st
taps1_coef = taps1_coef(1:14);
taps2_coef = taps2_coef(1:14);

% Start values for Register 2
StartValues = e5bq_init2(5*(PRN-1)+1:5*PRN);
StartValues = dec2bin(base2dec(StartValues,8))-48;

% Consider when the first bit in the start values for register 2 is 0:
% dec2bin() and base2dec() omit the higher bit of 0
Register2 = zeros(1,14);
Register2(end-length(StartValues)+1:end) = StartValues;

%% Generating the E5bQ code =========================================
Pri_E5bQ  = zeros(1,10230);
for ind = 1:length(Pri_E5bQ)
    % Register output mutiplied by the coefficeients
    RegOut1 = Register1 .* taps1_coef;
    RegOut2 = Register2 .* taps2_coef;
    
    % E5bQ primary code sequence output (bipolar format)
    Pri_E5bQ(ind) = (1 - 2* RegOut1(1)) * (1 - 2* RegOut2(1));
    
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
% ----- E5bQ primary or tiered code generation ----------------------------
if flag ==1
    E5bQ = Pri_E5bQ;
elseif  flag ==2
    % Secondary code
    SecondaryCode = generateE5bQ_secondary(PRN);
    
    % -- Combine the promary and secondary codes to form the tiered code --
    E5bQ = zeros(1,10230*100);
    for ind = 1:100
        % 100ms E5bQ tiered code
        E5bQ( (ind-1)*10230+1 : ind*10230 ) = Pri_E5bQ * SecondaryCode(ind);
    end
end