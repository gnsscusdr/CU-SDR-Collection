function [eph, firstSubFrame,TOD] = NAVdecoding(I_P_InputBits,settings)

% findPreambles finds the first preamble occurrence in the bit stream of
% each channel. The preamble is verified by check of the spacing between
% preambles (6sec) and parity checking of the first two words in a
% subframe. At the same time function returns list of channels, that are in
% tracking state and with valid preambles in the nav data stream.
%
%[eph, firstSubFrame,TOD] = CNAVdecoding(I_P_InputBits)
%
%   Inputs:
%       I_P_InputBits   - output from the tracking function
%
%   Outputs:
%       firstSubframe   - Starting positions of the first message in the 
%                       input bit stream I_P_InputBits in each channel. 
%                       The position is CNAV bit(20ms before convolutional decoding) 
%                       count since start of tracking. Corresponding value will
%                       be set to inf if no valid preambles were detected in
%                       the channel.
%       TOD             - Time Of Week (TOW) of the first message(in seconds).
%                       Corresponding value will be set to inf if no valid preambles
%                       were detected in the channel.
%       eph             - SV ephemeris. 

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Jakob Almqvist, Nagaraj C. Shivaramaiah and Dennis M. Akos
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
% $Id: findPreambles.m,v 1.1.2.10 2006/08/14 11:38:22 dpl Exp $


%--- Initialize ephemeris structute  --------------------------------------
% This is in order to make sure variable 'eph' for each SV has a similar 
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
eph = eph_structure_init();

% Starting positions of the first message in the input bit stream
firstSubFrame = inf;

% TOD of the first message
TOD = inf;

%% Bit and frame synchronization ====================================
% Preamble search can be delayed to a later point in the tracking results
% to avoid noise due to tracking loop transients
searchStartOffset = 0;

%--- Generate the preamble pattern ----------------------------------------
preamble_bits = [1 1 1 1 1 -1 -1 -1 1 1 -1 1 1 1 -1 1 -1 1 -1 ...
                                             -1 -1 -1 1 -1 -1 1 -1 1 1 -1];

% "Upsample" the preamble - make 10 vales per one bit. The preamble must be
% found with precision of a sample.
preamble_ms = kron(preamble_bits, ones(1, 10));

% Correlate tracking output with preamble ================================
% Read output from tracking. It contains the navigation bits. The start
% of record is skiped here to avoid tracking loop transients.
bits = I_P_InputBits(1 + searchStartOffset : end);

% Now threshold the output and convert it to -1 and +1
bits(bits > 0)  =  1;
bits(bits <= 0) = -1;

% Correlate tracking output with the preamble
tlmXcorrResult = xcorr(bits, preamble_ms);

% Find all starting points off all preamble like patterns ================
clear index
clear index2

xcorrLength = (length(tlmXcorrResult) +  1) /2;

%--- Find at what index/ms the preambles start ------------------------
index = find(...
    abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 271)' + ...
    searchStartOffset + 300;

%% Analyze detected preamble like patterns ================================
for i = 1:size(index) % For each occurrence
    
    %--- Find distances in time between this occurrence and the rest of
    %preambles like patterns. If the distance is 6000 milliseconds (one
    %subframe), the do further verifications by validating the parities
    %of two GPS words
    
    index2 = index - index(i);
    
    if (~isempty(find(index2 == 2000, 1)))
        
        %Be sure we have enough data for the decoding. 15 string are
        %needed. Each string represent 2000 IP values except for the last
        %one, for which we need only 1700 IP values
        if index(i) > settings.msToProcess - (14*2000 + 1700)
            disp('    Not enough data to decode this SV (subframe decoded too late)')
            return
        end

        
        %=== Re-read bit values for preamble verification =============
        % Preamble occurrence is verified by checking the parity of
        % the first two words in the subframe. Now it is assumed that
        % bit boundaries are known. Therefore the bit values over 10ms
        % are combined to increase receiver performance for noisy
        % signals.In total 85 bits must be read: 77 navigation bits and
        % 8 Hamming code bits. The index is pointing at the start of
        % the first navigation bit.
        biBinaryBits = I_P_InputBits(index(i): ...
            index(i) + 1700 -1)';
        
        %--- Combine the 10 values of each bit ------------------------
        biBinaryBits = reshape(biBinaryBits, 10, (size(biBinaryBits, 1) / 10));
        biBinaryBits = sum(biBinaryBits);
        
        % Now threshold and make it 0 and 1
        biBinaryBits(biBinaryBits > 0)  = 1;
        biBinaryBits(biBinaryBits <= 0) = 0;
        
        %--- Convert from bi-binary to relative code-------------------
        oddBits = biBinaryBits(1:2:169);
        evenBits = biBinaryBits(2:2:170);
        relRevBits = (oddBits - evenBits +1 )./2;
        
        %--- Convert from relative code to data sequence and checking
        % bits
        revBits(1) = 0;
        revBits(2:85) =  xor(relRevBits(1:84),relRevBits(2:85));
        
        %--- Reverse to correct bit order -----------------------------
        bits = revBits((85):-1:1);
        
        %--- Check the parity of the TLM and HOW words ----------------
        if (dataVerification(bits) ~= 0)
            % Verification OK. Record the TM start position. Skip
            % the rest of the TM pattern checking for this channel
            % and process next channel.
            if (firstSubFrame == inf)
                firstSubFrame = index(i);
            end
            
            %--- Copy 15 strings long record from tracking output ------------------
            % Skip the time marks.
            clear navBiBinaryBitsSamples
            for stringIndex = 1:15
                navBiBinaryBitsSamples((stringIndex-1)*1700 + 1: stringIndex*1700,1) = ...
                    I_P_InputBits( firstSubFrame + (stringIndex-1)*2000 : ...
                    firstSubFrame + (stringIndex-1)*2000 + (170 * 10 ) -1)';
            end
            
            %--- Ephemeris decoding ---------------------------------------------
            [eph, TOD] = ephemeris(navBiBinaryBitsSamples);
            break;
        end % if verification is OK ...
        
    end % if (~isempty(find(index2 == 2000)))
end % for i = 1:size(index)

% Return if no valid preamble was detected
if firstSubFrame == inf
    disp('Could not find valid preambles in this channel! ');
    return
end
