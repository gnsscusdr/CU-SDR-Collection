function [eph, firstSubFrame] = NAVdecodingv2(I_P,settings)

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
    
    
    %%  Initialize ephemeris structute  ---------------------------------------
    % Preamble search can be delayed to a later point in the tracking results
    % to avoid noise due to tracking loop transients
    searchStartOffset = 0;    
    % TOW of the first message
    TOW = inf;
    %Creates a cyclic redundancy code (CRC) detector System object
    % crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
    crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);
    
    % Convert CNAV-producing convolutional code polynomials to trellis description
    % Note that the difference from GPS is that the second branch G2 is
    % inverted at the end (see ICD)
    trellis = poly2trellis(7,[171 ~133]);
    % Viterbi traceback depth for vitdec(function)
    tblen = 35;
    % Secondary code is "842E9"
    % (Consider each bit from the nav message is a full secondary code.)
    secondCode = 1-2*[1 0 0 0 0 1 0 0 0 0 1 0 1 1 1 0 1 0 0 1];
    % Sync pattern 
    anti_sync_bits = [-1 1 -1 -1 1 -1 -1 -1 1 1 1 1];
    % Define eph structure and initialize validity flag to zero.
    eph = eph_structure_init();
    firstSubFrame = 0;
    eph.flag = 0;
    TOWflg = 1;
    wordCount = [];
    %% Find preamble to decode page -------------------------------------------
    % "Upsample" the preamble - make 20 vales per one bit. The preamble must be
    % found with precision of a sample.
    % This does: 
    % [ [1 0 1 1 0 1 1 1 0 0 0 0]*-1;
    %   [1 0 1 1 0 1 1 1 0 0 0 0]* 1  ];
    preamble_ms = kron(anti_sync_bits, secondCode);
    
    % Use the prompt correlator as the symbol stream, skip start of record if
    %   set in initSettings to avoid tracking loop transients
    symbols = I_P(1 + searchStartOffset : end);
    
    % Now threshold the output and convert it to -1 and +1
    symbols(symbols > 0)  =  1;
    symbols(symbols <= 0) = -1;
    
    % Correlate tracking output with the preamble
    tlmXcorrResult = xcorr(symbols, preamble_ms);
    xcorrLength = (length(tlmXcorrResult) +  1) /2;
    
        
    % Find at what index/ms the preambles start
    index = find(round(abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1))) >= 239.99)';
    % Check each index is 10e3 a part from at least another index
    newIndex = [];
    for idx=1:size(index,1)
        if find(abs(index-index(idx))==10e3)
            newIndex = [newIndex;index(idx)];
        end
    end
    index = newIndex;
    if isempty(index)
        fprintf('    Sync pattern not found in navigation data')
        return
    end
    % Save the starting point
    firstSubFrame = index(1);
    %% Start decoding current subframe, page by page. -------------------------
    %Indicator the requisite messages are all decoded
    wordType = 0;
    validCounter = 0;
    % Analyze detected preamble like patterns ================================
    for ind = 1:size(index) % For each occurrence
        % Check there are enough bits.
        if (settings.msToProcess - index(ind)) < 10e3 
            fprintf('    Not enough symbols to demodulate navData!', ind)
            break
        end
    
        % Select a whole page from the current subframe (1000 symbols, or 500 nav bits)
        page = symbols(index(ind):(index(ind) + 10e3 - 1));
        % Convert from symbols to navigation bits 
        navBits = (secondCode*reshape(page, 20, 500)) > 0;
        %--- Correct polarity of the all data bits according to preamble bits
        sync_bits = [1 0 1 1 0 1 1 1 0 0 0 0];
        if navBits(1:12)*sync_bits'~=6
            navBits = not(navBits);
        end
        % Remove preamble from data
        navBits = navBits(13:end);
    
        %--- De-interleave page parts -----------------------------------------
        symMat = reshape(navBits,61,8)';
        pageSym = reshape(symMat,1,[])';
        
        %--- Remove convolutional encoding from page parts --------------------
        decBits = vitdec(pageSym,trellis,tblen,'trunc','hard');
        
        %--- Retrieve full page --------------------------------------------
        page = decBits(1:238)';
        
        %--- Check the CRC ----------------------------------------------------
        [~,frmError] = step(crcDet,page');
        if (frmError)
%--- Decode the message type ------------------------------------------
            fprintf('\nChecksum error.\n')
            return 
        end
        
        %% Retrieve ephemeris data from current page --------------------------
        navWord = dec2bin(page)';
        [eph, validWord, wordType]= ephemeris(navWord, eph);
        % Correct TOW if could be retrieved
        if TOWflg && ~isempty(eph.TOW)
            % Correct TOW to time for first page part
            TOW = eph.TOW - (ind-1)*10;
            TOWflg = 0;
        end
        % Increase counter only if a different page has been read
        if ~any(wordType == wordCount)
            wordCount = [wordCount, wordType];
            validCounter = validCounter + validWord;
        end
        % Check if enough navigation pages have been decoded
        if validCounter ==4
            break
        end
    end %
    eph.TOW = TOW;
    eph.flag = floor(validCounter/4);        
end
