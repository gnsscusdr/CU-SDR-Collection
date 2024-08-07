function acqResults = acquisitionJMBF(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for Galileo E5b SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
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
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Condition input signal to speed up acquisition ===========================
% If input IF signal freq. is too high, a resampling strategy is applied
% to speed up the acquisition. This is user selectable.
if (settings.samplingFreq > settings.resamplingThreshold && ...
        settings.resamplingflag == 1)

    %--- Filiter out signal power outside the main lobe of CM code ------------------
    fs = settings.samplingFreq;
    IF = settings.IF;
    % Bandwidth of E5b mian lobe
    BW = 20.46e6;
    % Filter parameter
    w1 = (IF)-BW/2;
    w2 = (IF)+BW/2;
    wp = [w1*2/fs w2*2/fs];
    % Filter coefficients
    b  = fir1(700,wp);
    % Filter operation
    longSignal = filtfilt(b,1,longSignal);

    % --- Find resample frequency -----------------------------------------
    % Refer to bandpass sampling theorem(Yi-Ran Sun,Generalized Bandpass
    % Sampling Receivers for Software Defined Radio)

    % Upper boundary frequency of the bandpass IF signal
    fu = settings.IF + BW/2;

    % Lower freq. of the acceptable sampling Freq. range
    n = floor(fu/BW);
    if (n<1)
        n = 1;
    end
    lowerFreq = 2*fu/n;

    % Lower boundary frequency of the bandpass IF signal
    fl = settings.IF - BW/2;

    % Upper boundary frequency of the acceptable sampling Freq. range
    if(n>1)
        upperFreq = 2*fl/(n-1);
    else
        upperFreq = lowerFreq;
    end

    % Save orignal Freq. for later use
    oldFreq = settings.samplingFreq;

    % Take the center of the acceptable sampling Freq. range as
    % resampling frequency. As settings are used to generate local
    % CM code samples, so assign the resampling freq. to settings.
    % This can not change the settings.samplingFreq outside this
    % acquisition function.
    settings.samplingFreq = ceil((lowerFreq + upperFreq)/2);

    %--- Downsample input IF signal -------------------------------------------------
    % Signal length after resampling
    signalLen = floor((length(longSignal)-1) /oldFreq * settings.samplingFreq);
    % Resampled signal index
    index = ceil((0:signalLen-1)/settings.samplingFreq * oldFreq);
    index(1) = 1;
    % Resampled signal
    longSignal = longSignal(index);

    % For later use
    oldIF = settings.IF;

    % Resampling is equivalent to down-converting the original IF by integer
    % times of resampling freq.. So the IF after resampling is equivalent to:
    settings.IF = rem(settings.IF,settings.samplingFreq);

end % resampling input IF signals

%% Acquisition initialization =======================================

%--- Varaibles for coarse acquisition -------------------------------------
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of 2ms local carrier wave (1ms for local duplicate,
% the other 1ms for zero padding)
phasePoints = (0 : (samplesPerCode * 2 -1)) * 2 * pi * ts;

% Number of the frequency bins for the specified search band
numberOfFreqBins = round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1;
% Carrier frequency bins to be searched
coarseFreqBin = zeros(1, numberOfFreqBins);

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 50);
% E5bI code phases of detected signals
acqResults.codePhase    = zeros(1, 50);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 50);

%--- Varaibles for fine acquisition ---------------------------------------
% Number of the frequency bins for fine acquisition: use 5Hz step
NumOfFineBins = round(settings.acqSearchStep / 5) + 1;
% Carrier frequencies of the frequency bins
FineFreqBins     = zeros(1, NumOfFineBins);
% Search results of all frequency bins
FineResult = zeros(1,NumOfFineBins);
% Coherent integration for each code
sumPerCode = zeros(1,100);
% Find phase points of the local carrier wave
finePhasePoints = (0 : (100*samplesPerCode-1)) * 2 * pi * ts;

%--- Input signal power for GLRT statistic calculation --------------------
sigPower = sqrt(var(longSignal(1:samplesPerCode)) * samplesPerCode);

% Perform search for all listed PRN numbers ...
fprintf('(');
for PRN = settings.acqSatelliteList

    %% Coarse acquisition ===========================================
    % Generate all E5bI+E5bQ primary codes and sample them according to the sampling freq.
    E5bICodesTable = makeE5bITable(settings,PRN);
    E5bQCodesTable = makeE5bQTable(settings,PRN);

    % generate local code duplicate to do correlate
    localE5bICode = [E5bICodesTable, zeros(1,samplesPerCode)];
    localE5bQCode = [E5bQCodesTable, zeros(1,samplesPerCode)];

    % Search results of all frequency bins and code shifts (for one satellite)
    results = zeros(numberOfFreqBins, samplesPerCode*2);

    %--- Perform DFT of PRN code ------------------------------------------
    E5bICodeFreqDom = conj(fft(localE5bICode));
    E5bQCodeFreqDom = conj(fft(localE5bQCode));

    %--- Make the correlation for whole frequency band (for all freq. bins)
    for freqBinIndex = 1:numberOfFreqBins

        %--- Generate carrier wave frequency grid  -----------------------
        coarseFreqBin(freqBinIndex) = settings.IF + settings.acqSearchBand - ...
            settings.acqSearchStep * (freqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(-1i * coarseFreqBin(freqBinIndex) * phasePoints);

        %--- Do correlation -----------------------------------------------
        for nonCohIndex = 1: settings.acqNonCohTime
            % Take 2ms vectors of input data to do correlation
            signal = longSignal((nonCohIndex - 1) * samplesPerCode + ...
                1 : (nonCohIndex + 1) * samplesPerCode);
            % "Remove carrier" from the signal
            I      = real(sigCarr .* signal);
            Q      = imag(sigCarr .* signal);

            %--- Convert the baseband signal to frequency domain --------------
            IQfreqDom = fft(I + 1i*Q);

            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            convE5bI = IQfreqDom .* E5bICodeFreqDom;
            convE5bQ = IQfreqDom .* E5bQCodeFreqDom;

            %--- Perform inverse DFT and store correlation results ------------
            cohRresult = abs(ifft(convE5bI)) + abs(ifft(convE5bQ));
            % Non-coherent integration
            results(freqBinIndex, :) = results(freqBinIndex, :) + cohRresult;
        end % nonCohIndex = 1: settings.acqNonCohTime
    end % freqBinIndex = 1:numberOfFreqBins

    %% Look for correlation peaks for coarse acquisition ============
    % Find the correlation peak and the carrier frequency
    [~, acqCoarseBin] = max(max(results, [], 2));
    % Find code phase of the same correlation peak
    [peakSize, codePhase] = max(max(results));
    % Store GLRT statistic
    acqResults.peakMetric(PRN) = peakSize/sigPower/settings.acqNonCohTime;

    %% Fine resolution frequency search =============================
    % If the result is above threshold, then there is a signal ...
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        %--- Initialize arrays to speed up the code -----------------------
        acqResults.carrFreq(PRN) = coarseFreqBin(acqCoarseBin);
        % Save code phase acquisition result
        acqResults.codePhase(PRN) = codePhase;


        %% Downsampling recovery ============================================
        % Find acquisition results corresponding to orignal sampling freq
        if (exist('oldFreq', 'var') && settings.resamplingflag == 1)
            % Find code phase
            acqResults.codePhase(PRN) = floor((codePhase - 1)/ ...
                settings.samplingFreq * oldFreq)+1;

            % Doppler frequency
            doppler = acqResults.carrFreq(PRN) - settings.IF;

            % Carrier freq. corresponding to orignal sampling freq
            acqResults.carrFreq(PRN) = doppler + oldIF;
        end
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold

end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');