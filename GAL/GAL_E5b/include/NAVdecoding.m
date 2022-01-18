function [eph, firstSubFrame,TOW] = NAVdecoding(I_P)

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
%       firstSubFrame   - Starting positions of the first message in the
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
% $Id: findPreambles.m,v 1.1.2.10 2017/01/19 21:13:22 dpl Exp $


%--- Initialize ephemeris structute  --------------------------------------
% This is in order to make sure variable 'eph' for each SV has a similar
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
eph = eph_structure_init();

% Preamble search can be delayed to a later point in the tracking results
% to avoid noise due to tracking loop transients
searchStartOffset = 0;

% Starting positions of the first message in the input bit stream
firstSubFrame = inf;

% TOW of the first message
TOW = -1;

%Creates a cyclic redundancy code (CRC) detector System object
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Convert CNAV-producing convolutional code polynomials to trellis description
% Note that the difference from GPS is that the second branch G2 is
% inverted at the end (see ICD)
trellis = poly2trellis(7,[171 ~133]);

% Viterbi traceback depth for vitdec(function)
tblen = 35;

%--- Generate the sync pattern --------------------------------------------
% Secondary code is "E"
secondCode = [-1 -1 -1 1];

% The preamble is [0 1 0 1 1 0 0 0 0 0], and the antipodal form of the
% preamble pattern is:
preamble_bits = [1 -1 1 -1 -1 1 1 1 1 1];
sync_bits = preamble_bits <0;

% "Upsample" the preamble - make 4 vales per one bit. The preamble must be
% found with precision of a sample.
preamble_ms = kron(preamble_bits, secondCode);

% Use the prompt correlator as the symbol stream, skip start of record if
%   set in initSettings to avoid tracking loop transients
bits = I_P(1 + searchStartOffset : end);

% Now threshold the output and convert it to -1 and +1
bits(bits > 0)  =  1;
bits(bits <= 0) = -1;

% Correlate tracking output with the preamble
tlmXcorrResult = xcorr(bits, preamble_ms);

% Find all starting points off all preamble like patterns -----------------
clear index
xcorrLength = (length(tlmXcorrResult) +  1) /2;

% Find at what index/ms the preambles start
index = find(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 39.99)';
if isempty(index)
    eph = [];
    return
else
    % Analyze detected preamble like patterns ================================
    for ind = 1:size(index) % For each occurrence
        
        % Check distance to all other possible sync patterns for current start
        index2 = index - index(ind);
        %     index(ind)
        % Check spacing, need two 250 symbol pages (even and odd)
        if (~isempty(find(index2 == 250*4,1)))&& ...
                ( length(bits(index(ind):end))> 7500*4 )
            %=== Read bit vales for CRC-24Q check and ephemeris decoding ======
            % Search every possible preamble pattern.
            temp_bits = bits(index(ind):index(ind)+ 7500*4 -1);
            
            % Group every 10 I_P to a row for corresponding bit
            I_P_group = reshape(temp_bits,4,[])';
            
            % Wipe off the 2nd code and form data bits
            navBits = sum(I_P_group .* repmat(secondCode,7500,1),2)';
            navBits = navBits <0;
            %--- Correct polarity of the all data bits according to preamble bits
            if(~isequal(navBits(1:10),sync_bits))
                navBits = not(navBits);
            end
            
            % Decoding the even page pert -------------------------------------
            % Pull out implied pages from the detected preamble
            pageSymInt1 = navBits(11:250);
            
            % De-interleave symbols
            symMat1 = reshape(pageSymInt1,30,8)';
            pageSym1 = reshape(symMat1,1,[])';
            
            % Remove convolutional encoding from implied pages
            decBits1 = vitdec(pageSym1,trellis,tblen,'trunc','hard');
            
            % Remove tail bits
            decBits1 = decBits1(1:114);
            
            % The even page part is always the first one
            if (decBits1(1) ~= 0)
                continue
            end
            
            % Decoding the odd page pert --------------------------------------
            % Pull out implied pages from the detected preamble
            pageSymInt2 = navBits(261:500);
            
            % De-interleave symbols
            symMat2 = reshape(pageSymInt2,30,8)';
            pageSym2 = reshape(symMat2,1,[])';
            % Remove convolutional encoding from implied pages
            decBits2 = vitdec(pageSym2,trellis,tblen,'trunc','hard');
            
            % Remove tail bits and reserved 2 bits
            decBits2 = decBits2(1:106);
            
            % Page Type (the second bit) equal to 0 to indicate the nominal
            % page type
            if (decBits1(1) == 0) && (decBits1(2) == 0) && (decBits2(2) == 0)
                pageBits = [decBits1' decBits2'];
            else
                continue
            end
            
            % Detect errors in input data using CRC ---------------------------
            [~,frmError] = step(crcDet,pageBits');
            
            % CRC-24Q check was OK. Then to decode ephemeris message by message.
            if (~frmError)
                %--- Just save for first message ------------------------------
                % firstSubFrame is the starting positions of the first message
                % in the input bit stream dataBits
                firstSubFrame = index(ind) + searchStartOffset;
                
                %--- Ephemeris decoding ---------------------------------------
                eph = ephemeris(navBits',eph);
                if eph.flag
                    TOW = eph.TOW;
                end
                break
            end % if CRC is OK ...
            
        end % Check spacing, need two 250 symbol pages (even and odd)
    end
end % for ind = 1:size(index)
