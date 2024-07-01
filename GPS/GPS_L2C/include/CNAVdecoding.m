function [eph, firstSubFrame,TOW] = CNAVdecoding(I_P_InputBits)

% findPreambles finds the first preamble occurrence in the bit stream of
% each channel. The preamble is verified by check of the spacing between
% preambles (6sec) and parity checking of the first two words in a
% subframe. At the same time function returns list of channels, that are in
% tracking state and with valid preambles in the nav data stream.
%
%[eph, firstSubFrame,TOW] = CNAVdecoding(I_P_InputBits)
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
%       TOW             - Time Of Week (TOW) of the first message(in seconds).
%                       Corresponding value will be set to inf if no valid preambles
%                       were detected in the channel.
%       eph             - SV ephemeris. 

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

% TOW of the first message
TOW = inf;

%--- Take even number of input bits to do decoing ------------------------ 
evenLen = length(I_P_InputBits) - rem(length(I_P_InputBits),2);
dataBits = I_P_InputBits(1:evenLen);

% Now threshold the output and convert it to 0 and 1 
dataBits = (dataBits < 0);

% Convert CNAV-producing convolutional code polynomials to trellis description
trellis = poly2trellis(7,[171 133]);

% Viterbi traceback depth for vitdec(function)
tblen = 35;

%--- Generate the preamble pattern ----------------------------------------
preamble_01   = [1 0 0 0 1 0 1 1];
% Antipodal form of the preamble pattern
preamble_bits = [-1 1 1 1 -1 1 -1 -1];

%% CNAV data decoding =====================================================
for G1orG2 = 1:2
    % The first bit in dataBits may be G1 or G2 ouput, the bit
    % stream to be decoded must start at bit of G1 output. So we must
    % search the first two bits to find the right one corresponding to G1.
    decodedBits = vitdec(dataBits(G1orG2:end-(G1orG2-1)),...
                                             trellis,tblen,'trunc','hard');

    % Now threshold the output and convert it to -1 and +1
    antipodalBits = 1 - 2*decodedBits;
    
    % Correlate tracking output with the preamble
    tlmXcorrResult = xcorr(antipodalBits, preamble_bits);
    
    %% Find all starting points off all preamble like patterns ==============
    clear index
    
    xcorrLength = (length(tlmXcorrResult) +  1) /2;
    
    %--- Find at what index/ms the preambles start ------------------------
    index = find(...
        round(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1))) == 8.00)';
    % Analyze detected preamble like patterns ================================
    for i = 1:height(index) % For each occurrence
        
        % Ensure the search for i has the number of a whole message(300) 
        if ((length(decodedBits) - index(i) +1) >=300 )
            
            %=== Read bit vales for CRC-24Q check and ephemeris decoding ============
            % Search every possible preamble pattern. This can allow decoding 
            % every message even when the data bit stream has only one message
            % without bit error.
            navBits = decodedBits(index(i):end);
            
            %--- Correct polarity of the all data bits 
            % according to preamble bits --------------------------------------------
            if(~isequal(navBits(1:8),preamble_01))
                navBits = not(navBits);
            end
            
            %--- To do CRC-24Q check of current mesage ------------------------------
            %Creates a cyclic redundancy code (CRC) detector System object
            % The generation polynomial coefficients
            crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
            
            % Detect errors in input data using CRC
            [~,frmError] = step(crcDet,navBits(1:300)');          
            
            % CRC-24Q check was OK. Then to decode ephemeris message by message.
            if (~frmError) 
                %--- Convert from decimal to binary -----------------------------------
                % The function ephemeris_L2CM expects input in binary form. In Matlab 
                % it is a string array containing only "0" and "1" characters.
                navBitsBin = dec2bin(navBits);
                %--- Ephemeris decoding ---------------------------------------------
                [eph, TOW_temp] = ephemeris(navBitsBin',eph);
                %--- Just save for first message ------------------------------------
                % firstSubFrame is the starting positions of the first message
                % in the input bit stream dataBits
                if (firstSubFrame == inf)
                    firstSubFrame = index(i) * 2 + G1orG2 - 2;
                    TOW = TOW_temp;
                    eph.TOW = TOW;
                end
            end % if CRC is OK ...
            
        end % if (~isempty(find(index2 == 300)))
    end % for i = 1:size(index)
    
    % If the first bit in dataBits is just G1 output, then quit the loop
    if (firstSubFrame ~= inf)
        break
    end
end % for G1orG2 = 1:2
