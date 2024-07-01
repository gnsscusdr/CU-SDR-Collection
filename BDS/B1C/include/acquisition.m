function acqResults = acquisition(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 20 ms of raw IF signal from the front-end..
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
% (C) Developed for BDS B1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%
% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and
% implementation of an open-source BDS-3 B1C/B2a SDR receiver.
% GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
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


%% Condition input signal to speed up acquisition ===================

% If input IF signal freq. is too high, a resampling strategy is applied
% to speed up the acquisition. This is user selectable.
if (settings.samplingFreq > settings.resamplingThreshold && ...
        settings.resamplingflag == 1)

    %--- Filiter out signal power outside the main lobe of CM code --------
    fs = settings.samplingFreq;
    IF = settings.IF;
    % Bandwidth of CM mian lobe
    BW = 9e6;
    % Filter parameter
    w1 = (IF)-BW/2;
    w2 = (IF)+BW/2;
    wp = [w1*2/fs-0.002 w2*2/fs+0.002];
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

    %--- Downsample input IF signal ---------------------------------------
    % Signal length after resampling
    signalLen = floor((length(longSignal)-1) /oldFreq * settings.samplingFreq);
    % Resampled signal index
    index = ceil((0:signalLen-1)/settings.samplingFreq * oldFreq);
    index(1) = 1;
    % Resampled signal
    longSignal = longSignal(index);

    % Foe latter use
    oldIF = settings.IF;

    % Resampling is equivalent to down-converting the original IF by integer
    % times of resampling freq.. So the IF after resampling is equivalent to:
    settings.IF = rem(settings.IF,settings.samplingFreq);

end % resampling input IF signals

%% Acquisition initialization =======================================
%--- Varaibles for coarse acquisition -------------------------------------
% Number of samples per B1C code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Number of samples of settings.acqCohT ms B1C spreading codes
samplesXmsLen= round(samplesPerCode /10*settings.acqCohT);
% Number of samples of 10 plus acqCohT(X) ms B1C spreading codes used to do
% zero padding correlation for coarse acquisition
len10PlusXms = round(samplesPerCode /10*(10+settings.acqCohT));
% Cut 10 plus X cm input signal to do correlate
sig10PlusXms = longSignal(1:len10PlusXms);
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of the local carrier wave
phasePoints = (0 : (len10PlusXms-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band
numberOfFrqBins = round(settings.acqSearchBand * 2 / settings.acqStep) + 1;
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, len10PlusXms);
% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);

%--- Initialize acqResults and related variables --------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, max(settings.acqSatelliteList));
% PRN code phases of detected signals
acqResults.codePhase    = zeros(1, max(settings.acqSatelliteList));
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, max(settings.acqSatelliteList));

%--- Varaibles for fine acquisition ---------------------------------------
% Number of the frequency bins for fine acquisition: use
% settings.acqStep*2 fine acquisition band, and 25Hz step
fineStep = 25;
NumOfFineBins = round(settings.acqStep/25)*2 + 1;
% Carrier frequencies of the frequency bins
FineFrqBins     = zeros(1, NumOfFineBins);
% Non-coherent results of all frequency bins
FineResult = zeros(1,NumOfFineBins);
% Phase points of the local carrier wave
finePhasePoints = (0 : samplesPerCode-1) * 2 * pi * ts;

%--- Input signal power for GLRT statistic calculation --------------------
sigPower = sqrt(var(sig10PlusXms(1:samplesXmsLen)) * samplesXmsLen);
%--- Generate carrier wave frequency grid  -----------------------
initFreq = settings.IF + settings.acqSearchBand;
% Generate local sine and cosine
sigCarr = exp(-1i*initFreq* phasePoints);
% "Remove carrier" from the signal
I1      = real(sigCarr .* sig10PlusXms);
Q1      = imag(sigCarr .* sig10PlusXms);
% Convert the baseband signal to frequency domain
IQfreqDom = fft(I1 + 1i*Q1);
% Perform search for all listed PRN numbers ...
fprintf('(');
for PRN = settings.acqSatelliteList
    %% Coarse acquisition ===========================================

    % Generate B1C data codes and sample them according to the sampling freq.
    DataPriTable = makeDataTable(settings,PRN);
    % generate local code duplicate to do correlate
    localData = [DataPriTable(1:samplesXmsLen) ...
        zeros(1,len10PlusXms - samplesXmsLen)];

    % Perform DFT of B1C data code
    DataPriFreqDom = conj(fft(localData));

    % Pilot signal acquisition
    if (settings.pilotACQflag == 1)
        PilotPriTable = makePilotTable(settings,PRN);
        localPilot = [PilotPriTable(1:samplesXmsLen) ...
            zeros(1,len10PlusXms - samplesXmsLen)];
        PilotPriFreqDom = conj(fft(localPilot));
    end
    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins
        IQfreqDomShifted = circshift(IQfreqDom, frqBinIndex - 1);
        % Multiplication in the frequency domain (correlation in time domain)
        convCodeIQ1 = IQfreqDomShifted .* DataPriFreqDom;
        % Perform inverse DFT and store correlation results
        results(frqBinIndex, :) = abs(ifft(convCodeIQ1));

        % Pilot signal acquisition
        if (settings.pilotACQflag == 1)
            convCodeIQ1 = IQfreqDomShifted .* PilotPriFreqDom;
            % Non-coherent combining of data and pilot results
            results(frqBinIndex, :) = (results(frqBinIndex, :)* sqrt(11)+ ...
                abs(ifft(convCodeIQ1))* sqrt(29) )/ sqrt(40);
        end

    end % frqBinIndex = 1:numberOfFrqBins

    %% Look for correlation peaks in the results ====================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak

    % Find the correlation peak and the carrier frequency
    [~, frequencyBinIndex] = max(max(results, [], 2));
    selFreq = initFreq - (frequencyBinIndex -1)*settings.acqStep;
    % Find code phase of the same correlation peak
    [peakSize, codePhase] = max(max(results));
    % Store GLRT statistic
    acqResults.peakMetric(PRN) = peakSize/sigPower;

    % To prevent index from exceeding matrix dimensions in the fine
    % qacquisition, move to previous code start position.
    if (codePhase+samplesPerCode-1) > length(longSignal)
        codePhase = codePhase - samplesPerCode;
    end

    % If the result is above threshold, then there is a signal ...
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        %% Fine resolution frequency search =========================
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        % Cut 10ms input signal
        signal0DC = longSignal(codePhase:(codePhase + samplesPerCode-1));
        % Remove B1C code modulation from the original signal
        xCarrier = signal0DC .* DataPriTable;

        % Pilot signal acquisition
        if (settings.pilotACQflag == 1)
            xCarrierPilot = signal0DC .* PilotPriTable;
        end

        %--- Search different frequency bins ------------------------------
        for FineBinIndex = 1 : NumOfFineBins
            % Carrier frequencies of the fine frequency bins
            FineFrqBins(FineBinIndex) = selFreq + settings.acqStep - fineStep * (FineBinIndex - 1);
            % Generate local sine and cosine
            sigCarr10cm = exp(-1i*FineFrqBins(FineBinIndex) * finePhasePoints);
            FineResult(FineBinIndex) = abs(sum(xCarrier .* sigCarr10cm));

            %--- Pilot signal acquisition
            if (settings.pilotACQflag == 1)
                FineResult(FineBinIndex) = (FineResult(FineBinIndex)*11 + ...
                    abs(sum(xCarrierPilot .* sigCarr10cm))*29)/40;
            end
        end

        % Find the fine carrier freq. -------------------------------------
        % Corresponding to the largest noncoherent power
        [~,maxFinBin] = max(FineResult);
        acqResults.carrFreq(PRN) = FineFrqBins(maxFinBin);

        %signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end
        acqResults.codePhase(PRN) = codePhase;

        %% Downsampling recovery  ===================================
        % Find acquisition results corresponding to orignal sampling freq
        if (exist('oldFreq', 'var') && settings.resamplingflag == 1)
            % Find code phase
            acqResults.codePhase(PRN) = floor((codePhase - 1)/ ...
                settings.samplingFreq * oldFreq)+1;

            % Doppler frequency
            if (settings.IF >= settings.samplingFreq/2)
                % In this condition, the FFT computed freq. is symmetric
                % with the true frequemcy about half of the sampling
                % frequency, so we have the following:
                IF_temp = settings.samplingFreq - settings.IF;
                doppler = IF_temp - acqResults.carrFreq(PRN);
            else
                doppler = acqResults.carrFreq(PRN) - settings.IF;
            end
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
