function [flag,decodedBits] = BCH21_6Decoding(bits)

% This function decodes the BCH(21,6) data symbols of the first subframe of
% in the BCNAV1 message. A total of 2^6-1 hypotheses are repeated to find
% the 6 original symbols.
%
% [flag,outputBits] = BCH21_6Decoding(bits)
%
%   Inputs:
%       bits         - Input row vector of message symbols to be decoded
%                      with 21 bits in bipolar (-1, +1)
%   Outputs:
%       flag         - Indicator for  decoding of the BCH(21,6) data 
%                      symbols, 1:success, 0: fail.
%       decodedBits  - Decoded message symboles for the BCH(21,6)

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Written by Yafeng Li
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
%$Id: generateE5BIcode.m,v 1.1.2.5 2018/07/24 22:00:00 dpl Exp $

% Threshold for determining the success or fail of the decoding. This value
% can be reconfigued according to sepecific requirements
threshold = 20;

% Total 2^6-1 hypotheses of the original symbols before BCH(21,6) encoding 
hypoBits = zeros(2^6,6);

% Symbols after BCH(21,6) encoding for each hypothesis
encoded = zeros(1,21);

% Correlation values between the recieved data bits and the symbols after
% BCH(21,6) encoding for each hypothesis
correValue = zeros(1,2^6);

% Feedback positions of the register for BCH(21,6) encoding
reg1_FeedbackPos = [2 4 5 6];

%% Decoding for each hypothesis =====================================
for  hypothesis = 1:2^6
    % Generation of each hypothesis symbols -------------------------------
    % Change current hypothesis into binary format of "1" and "0"
    register1 = dec2bin(hypothesis-1) - 48;
    register1 = [zeros(1,6-length(register1)) register1]; %#ok<AGROW>
    
    % Save the binary format of current hypothesis symbols
    hypoBits(hypothesis,:) = register1;
    
    % Convert symbol polarity: 0 -> 1 and 1 -> -1
    register1(register1 == 1) = -1;
    register1(register1 == 0) = 1;
    
    % Flip "register1" left to right
    register1 = fliplr(register1);

    % Encoding of the BCH(21,6) -------------------------------------------
    for ind = 1:21
        encoded(ind) = register1(end);
        % Exclusive OR operation for feedback
        feedback1 = prod(register1(reg1_FeedbackPos));
        % Shift the register to right by one element
        register1 = circshift(register1',1)';
        register1(1) = feedback1;
    end
    % Compute the correlation values
    correValue(hypothesis) = sum(encoded .* bits);
end

% Determine the results according to the correlation values ---------------
[maxValue,pos] = max(correValue);

if maxValue >= threshold
    flag = 1;
    decodedBits = hypoBits(pos,:);
else
    flag = 0;
    decodedBits = [];
end
