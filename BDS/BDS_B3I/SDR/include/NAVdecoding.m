function [eph, subFrameStart,SOW] = NAVdecoding(I_P_InputBits,PRN,settings)

% findPreambles finds the first preamble occurrence in the bit stream of
% each channel. The preamble is verified by check of the spacing between
% preambles (6sec) and parity checking of the first two words in a
% subframe. At the same time function returns list of channels, that are in
% tracking state and with valid preambles in the nav data stream.
%
%[eph, subFrameStart,SOW] = CNAVdecoding(I_P_InputBits)
%
%   Inputs:
%       I_P_InputBits   - output from the tracking function
%
%   Outputs:
%       subFrameStart   - Starting positions of the first message in the
%                       input bit stream I_P_InputBits in each channel.
%                       The position is CNAV bit(20ms before convolutional decoding)
%                       count since start of tracking. Corresponding value will
%                       be set to inf if no valid preambles were detected in
%                       the channel.
%       SOW             - Time Of Week (SOW) of the first message(in seconds).
%                       Corresponding value will be set to inf if no valid preambles
%                       were detected in the channel.
%       eph             - SV ephemeris.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B3I SDR by Daehee Won, Yafeng Li, 
% Nagaraj C. Shivaramaiah and Dennis M. Akos. 
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%--------------------------------------------------------------------------
%
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

% CVS record:
% $Id: findPreambles.m,v 1.1.2.10 2017/01/19 21:13:22 dpl Exp $


%--- Initialize ephemeris structute  --------------------------------------
% This is in order to make sure variable 'eph' for each SV has a similar
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
eph = eph_structure_init();

% Starting positions of the first message in the input bit stream
subFrameStart = inf;

% SOW of the first message
SOW = inf;

%% Bit and frame synchronization ====================================
% Preamble search can be delayed to a later point in the tracking results
% to avoid noise due to tracking loop transients
searchStartOffset = 1000;

%--- Generate the preamble pattern ----------------------------------------
preamble_bits = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
NH_code = [-1 -1 -1 -1 -1 1 -1 -1 1 1 -1 1 -1 1 -1 -1 1 1 1 -1];

% "Upsample" the preamble - make 20 vales per one bit (D1) & 2 values per
% one bit (D2). The preamble must be found with precision of a sample.
preamble_D1 = kron(preamble_bits, -NH_code);
preamble_D2 = kron(preamble_bits, ones(1, 2));

%% Correlate tracking output with preamble ================================
% Read output from tracking. It contains the navigation bits. The start
% of record is skiped here to avoid tracking loop transients.
bits = I_P_InputBits(1 + searchStartOffset : end);

% Now threshold the output and convert it to -1 and +1
bits(bits > 0)  =  1;
bits(bits <= 0) = -1;

if ((1 <= PRN) && (PRN <=5 )) || ((59 <= PRN) && (PRN <= 63))
    % GEO satellites. Not modulated with NM code.
    % Correlate tracking output with the preamble
    tlmXcorrResult = xcorr(bits, preamble_D2);
    codrPerD = 2;
elseif (6 <= PRN) && (PRN <= 58)
    % MEO/IGSO satellites. Not modulated with NM code.
    % Correlate tracking output with the preamble
    tlmXcorrResult = xcorr(bits, preamble_D1);
    codrPerD = 20;
end

% Find all starting points off all preamble like patterns ================
clear index
clear index2

xcorrLength = (length(tlmXcorrResult) +  1) /2;

%--- Find at what index/ms the preambles start ------------------------
index = find(...
    abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) >= codrPerD*10)' + ...
    searchStartOffset;

%--- Making sure we have enough data for decoding
index = index(index < settings.msToProcess - 1500*20 + 300*codrPerD);
 
% Analyze detected preamble like patterns =============================
for i = 1:size(index) % For each occurrence
    
    %--- Find distances in time between this occurrence and the rest of
    %preambles like patterns. If the distance is 6000 milliseconds (one
    %subframe), the do further verifications by validating the parities
    %of two GPS words
    
    index2 = index - index(i);
    
    if (~isempty(find(index2 == 300*codrPerD, 1)))
        
        %=== Re-read bit vales for preamble verification ==============
        % Preamble occurrence is verified by checking the parity of
        % the first two words in the subframe. Now it is assumed that
        % bit boundaries a known. Therefore the bit values over 20ms are
        % combined to increase receiver performance for noisy signals.
        % in Total 62 bits mast be read :
        % 2 bits from previous subframe are needed for parity checking;
        % 60 bits for the first two 30bit words (TLM and HOW words).
        % The index is pointing at the start of TLM word.
        bits = I_P_InputBits(index(i): index(i) + 30 * codrPerD -1)';
        
        %--- Combine the 20 values of each bit ------------------------
        bits = reshape(bits, codrPerD, (size(bits, 1) / codrPerD));
        bits = sum(bits);
        
        % Now threshold and make it -1 and +1
        bits(bits > 0)  = 1;
        bits(bits <= 0) = 0;
        
        % Create a Galois field array from bits 
        bits = gf(bits,1);
        bits = bits(16:30);
        
        % Decode the received signal in code using an [15,11] BCH decoder
        [~,cnumerr] = bchdec(bits,15,11);
        
        %--- Check the parity of the TLM and HOW words ----------------
        if cnumerr == 0
            % Parity was OK. Record the preamble start position. Skip
            % the rest of preamble pattern checking for this channel
            % and process next channel.
            
            subFrameStart = index(i);
            break;
        end % if parity is OK ...
        
    end % if (~isempty(find(index2 == 6000)))
end % for i = 1:size(index)

% Exclude channel from the active channel list if no valid preamble was
% detected
if subFrameStart == inf
    disp(['    Could not find valid preambles for PRN # ',num2str(PRN),'!']);
    return
end

%% Decode ephemerides ===============================================
demod_bit = kron(ones(1,300*5), NH_code);

%=== Convert tracking output to navigation bits =======================
% GEO satellites. Not modulated with NM code.
if ((1 <= PRN) && (PRN <=5 )) || ((59 <= PRN) && (PRN <= 63))
    %--- Copy 5 sub-frames long record from tracking output ---------------
    % 1 frame has 1500 bits (5 subframes) and 1 bit is 2ms.
    navBitsSamples = I_P_InputBits(subFrameStart: ...
        subFrameStart + (1500 * 20) -1)';
    
    %--- Group every 2 vales of bits into columns ------------------------
    % 1 bit is 2ms long. Raw bit is converted to Navigation bit.
    % 3,000 raw bits --> 1500 Nav bits
    navBitsSamples = reshape(navBitsSamples, ...
        2, (size(navBitsSamples, 1) / 2));
    
    % MEO/IGSO satellites. Not modulated with NM code.
elseif (6 <= PRN) && (PRN <= 58)
    %--- Copy 5 sub-frames long record from tracking output ---------------
    % 1 frame has 1500 bits (5 subframes) and 1 bit is 20ms.
    navBitsSamples = I_P_InputBits(subFrameStart: ...
        subFrameStart + (1500 * 20) -1)';
    
    % Demodulating NH code
    navBitsSamples = navBitsSamples' .* demod_bit;
    
    %--- Group every 20 vales of bits into columns ------------------------
    % 1 bit is 20ms long. Raw bit is converted to Navigation bit.
    % 30,000 raw bits --> 1500 Nav bits
    navBitsSamples = reshape(navBitsSamples, ...
        20, (size(navBitsSamples, 2) / 20));
end

%--- Sum all samples in the bits to get the best estimate -------------
navBits = sum(navBitsSamples);

%--- Now threshold and make 1 and 0 -----------------------------------
% The expression (navBits > 0) returns an array with elements set to 1
% if the condition is met and set to 0 if it is not met.
navBits = (navBits > 0);

%--- Convert from decimal to binary -----------------------------------
% The function ephemeris expects input in binary form. In Matlab it is
% a string array containing only "0" and "1" characters.
%     navBitsBin = dec2bin(navBits);

%--- Decode data and extract ephemeris information ---
% 'ephemeris' is different from GPS.
[eph, SOW] = ephemeris(navBits, PRN);
