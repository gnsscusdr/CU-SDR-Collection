function B3Icode = generateB3Icode(PRN)
% generateB3Icode.m generates one of the Beidou satellite B3I codes.
%
% B3Icode = generateB3Icode(PRN)
%
%   Inputs:
%       PRN         - PRN number from #1 to #63 of the sequence. 
%       settings    - receiver settings
%
%   Outputs:
%       B3Icode      - a vector containing the desired B3I code sequence 
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
%$Id: generateB3Icode.m,v 1.1.2.5 2020/01/16 23:00:00 dpl Exp $


%--- B3I code length ------------------------------------------------------
CodeLength = 10230;

%--- Generate CA codes ----------------------------------------------------
% Feedback position of exclusive OR operation for CA codes
ca_FeedbackPos = [1 3 4 13];

% Initial state of CA register: 1 for 0, and -1 for 1, to perform exclusive 
% OR operation by multiply. 
ca_reg = ones(1,13) * -1;

% CA codes output
CA = zeros(1,CodeLength);

% CA register will be reset to all 1 (here -1 represents 1)
reset_state = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1];

for ind = 1:CodeLength
    CA(ind) = ca_reg(end);
    if isequal(ca_reg,reset_state)
        ca_reg = ones(1,13) * -1;
    else
        % Exclusive OR operation for feedback
        feedback = prod(ca_reg(ca_FeedbackPos));
        % shift the register to right by one element
        ca_reg = circshift(ca_reg',1)';
        ca_reg(1) = feedback;        
    end
end

%--- Generate CB codes ----------------------------------------------------
% Initial-state values from Beidou B3I ICD from PRN #1 to #63
B3I_init = [4;    11;   13;   22;   30;   36;   44;   48;   88;   104; ...
            116;  129;  376;  418;  458;  682;  696;  707;  1078; 2069; ...
            2248; 2574; 2596; 2731; 4294; 4436; 4647; 4978; 4986; 1; ...
            5209; 5539; 6061; 6488; 7130; 7165; 7403; 5879; 1681; 5080; ... 
            5938; 3983; 6208; 7223; 2996; 1814; 6906; 6144; 4713; 7406; ... 
            7264; 1766; 5347; 3515; 7951; 7054; 3884; 6067; 4230; 3803; ... 
            869; 3683; 1205];
        
% Feedback position of exclusive OR operation for CB codes
cb_FeedbackPos = [1 5 6 7 9 10 12 13];

% Initial state of CB register: 1 for 0, and -1 for 1, to perform exclusive 
% OR operation by multiply
cb_reg = ones(1,13) * -1;

% CB codes output
CB = zeros(1,CodeLength);

% CB code Advance(Chips) for the start of the CB code
resetPos = B3I_init(PRN);

% CB initial states
for ind = 1:resetPos
    % Exclusive OR operation for feedback
    feedback = prod(cb_reg(cb_FeedbackPos));
    % shift the register to right by one element
    cb_reg = circshift(cb_reg',1)';
    cb_reg(1) = feedback; 
end

% CB output
for ind = 1:CodeLength
    CB(ind) = cb_reg(end);
    % Exclusive OR operation for feedback
    feedback = prod(cb_reg(cb_FeedbackPos));
    % shift the register to right by one element
    cb_reg = circshift(cb_reg',1)';
    cb_reg(1) = feedback; 
end

%---  B3I PRN code output -------------------------------------------------
B3Icode = CB .* CA;