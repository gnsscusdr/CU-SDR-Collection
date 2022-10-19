function [eph, firstSubFrame,TOW] = NAVdecodingv2(I_P)

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
TOW = inf;

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
% Check each index is 1e3 a part from at least another index
newIndex = [];
for idx=1:size(index,1)
    if find(abs(index-index(idx))==1e3)
        newIndex = [newIndex;index(idx)];
    end
end
index = newIndex;
if isempty(index)
    fprintf('    Sync pattern not found in navigation data')
    return
end
subframe_sec = 15;
page_sym = 500;
validCount = 0;
TOWflg = 1;
wordlst = [];
%--- Reset search if IODNAV is different for any word
IODlst = [];
eph = eph_structure_init();
pagecnt = 0;
% Get Num of prompts required to decode a page
NsymSubFrameEnc = page_sym*length(secondCode)*1;
% Analyze detected preamble like patterns ================================
for ind = 1:size(index) % For each occurrence
        
        % Check there is enough data starting from the sync pattern
        if length(bits(index(ind):end)) > NsymSubFrameEnc
            %=== Read bit vales for CRC-24Q check and ephemeris decoding ======
            % Search every possible preamble pattern within the subframe
            % period
            page_bits = bits(index(ind):index(ind)+page_sym*length(secondCode) -1);
            
            % Because of 2dary code, each nav symbol is 4 bits
            I_P_group = reshape(page_bits,4,[])';
            
            % Wipe off the 2nd code and form data bits
            navBits = I_P_group*secondCode'<0;
            %--- Correct polarity of the all data bits according to preamble bits
            if sync_bits*navBits(1:10) < 3
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
            
            % Remove tail bits (1 even/odd + 1 page type + 120 data) 
            decBits1 = decBits1(1:114);
            
            % The even page part is always the first one
            if ~(decBits1(1) == 0)
                continue
            end
            
            % Decoding the odd page part --------------------------------------
            % Pull out implied pages from the detected preamble
            pageSymInt2 = navBits(261:end);
            
            % De-interleave symbols
            symMat2 = reshape(pageSymInt2,30,8)';
            pageSym2 = reshape(symMat2,1,[])';
            % Remove convolutional encoding from implied pages
            decBits2 = vitdec(pageSym2,trellis,tblen,'trunc','hard');
            
            % Remove tail bits and reserved 2 bits
            decBits2 = decBits2(1:106);
            
            % Page Type (the second bit) equal to 0 to indicate the nominal
            % page type
            if (decBits1(2) + decBits2(2)) == 0
                pageBits = [decBits1' decBits2'];
            else
                continue
            end
            
            % Detect errors in input data using CRC ---------------------------
            [~,frmError] = step(crcDet,pageBits');
            pagecnt = pagecnt + 1;
            % CRC-24Q check was OK. Then to decode ephemeris message by message.
            if (~frmError)
                %--- Just save for first message ------------------------------
                % firstSubFrame is the starting positions of the first message
                % in the input bit stream dataBits
                firstSubFrame = index(ind) + searchStartOffset;
                %---              
                navWordDec = [pageBits(3:114) pageBits(116:131)];
                %--- Convert navigation data word to binary ---------------------------
                navWord = dec2bin(navWordDec)';
                %--- Decode the message type ------------------------------------------
                wordType = bin2dec(navWord(1:6));
                %--- Check issue of data remains the same
                if wordType
                    [eph, validWord] = ephemerisv2(navWord, wordType, eph);
                    validCount = validCount + validWord;
                    %--- Add IODe to the list
                    if any([1:4] == wordType)
                        IODlst = [IODlst eph.IODnav];
                        %The moment a different IOD is received, reset count 
                        if (length(unique(IODlst))) > 1
                            validCount = 1;
                            IODlst = [eph.IODnav];
                            pagecnt = 1;
                        end
  
                    end
                    %--- Update TOW when read
                    if isempty(eph.TOW)
                        % Correct TOW to time for first page part 
                        TOW = eph.TOW - (pagecnt*2);
                        % Reset parameter in case it has to be re-read
                        eph.TOW = [];
                    end
                    %--- Check every 6 counts if IODs meet requirement
                    if validCount == 6
                        eph.flag = 1;
                        break                        
                    end   
                end % WordType is suitable for reading
             end % if CRC is OK ...
        else 
            % Not enough symbols to decode this page.
            disp('Not enough symbols to decode current page.')
        end % If enough symbols to decode page
end % For loop iterating sync patterns
   
end % for ind = 1:size(index)
