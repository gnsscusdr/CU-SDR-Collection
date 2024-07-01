function acqResults = acquisitionDMA(longSignal, settings)
%% Initialization ===ts======================================================
% Codes per block
Ncodes = 2;
% N blocks (Codes + zeros)
Nblocks = 4;
% Find number of samples per spreading code
samplesPerBlock = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / (Nblocks* settings.codeLength)));

% Load as many signal blocks as defined by the NonCoh sums
signal1= longSignal(1:samplesPerBlock);
signal2= longSignal(samplesPerBlock + 1: samplesPerBlock*2);
% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave
phasePoints = (0 : (samplesPerBlock-1)) * 2 * pi * ts;

% Find the current frequency resolution of the signal given the sampling frequency and a
% code length (BW/Nsamples)
freqResolution = settings.samplingFreq/samplesPerBlock;

% Find the number of frequency bins for the given frequency resolution and
% search band. 
numberOfFrqBins = round(settings.acqSearchBand*1e3/freqResolution) + 1;

% Measure or obtain the required frequency resolution for good correlation results 
% and see if it matches the one obtained above.

if isempty(settings.stepSize)  
    % If a stepSize is not given, it is calculated so as to capture half a
    % cycle of the signal for good correlation results.
    stepSize = 0.5/(Nblocks*settings.codeLength/settings.codeFreqBasis);
elseif settings.stepSize == freqResolution
    stepSize = settings.stepSize;
else
    % Define possible stepSizes that yield an integer number of sub-bins
    steps = [1:0.25:freqResolution/2];
    % Select just valid steps with zero remainder
    steps = steps(rem(freqResolution, [1:0.25:freqResolution/2]) == 0);
    % Chose the closest integer divisor
    stepDiff = steps-settings.stepSize;
    [~, minDiv] = min(abs(stepDiff));
    % Ensure the resulting bin size is equal or smaller (more restrictive)
    if stepDiff(minDiv) > 0
        stepSize = steps(minDiv -1);
    else
        stepSize = steps(minDiv);
    end
end

% Check if the required stepSize meets the resolution desired in settings or if
% interleaving is required
Nshifts = freqResolution/stepSize;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTableDMA(settings, Ncodes);

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 58);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 58);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 58);
% Initialize the starting frequency
initFreq = settings.IF + (settings.acqSearchBand/2) * 1000;
% Initialize matrix that accumulates correlation results
accMat = zeros((numberOfFrqBins - 1)*Nshifts +1,samplesPerBlock);

% Start acquisition 
fprintf('(');
idx = 0;
% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList
    tstart3 = tic;
    idx = idx + 1;
    codePhase = 0;
    freqShift = 0;
    prevmax = 0;
    corrVec = zeros(1, samplesPerBlock);
    AccCorrVec = zeros(1, samplesPerBlock);
    AccFreqVec = zeros(1, samplesPerBlock*Nshifts);
    %% Correlate signals ======================================================
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft([caCodesTable(PRN, :) zeros(1,samplesPerBlock/Ncodes)]));
    % Iterate over each sub-bin for better freq resolution
    for binIter = 1:Nshifts
        %--- Generate carrier wave frequency grid (1kHz step per 1msec) -----------
        initFreqShift = initFreq + (binIter-1)*(freqResolution/Nshifts);      
        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(-1i*initFreqShift*phasePoints);

            %--- "Remove carrier" from the signal -----------------------------
        I1      = real(sigCarr .* signal1);
        Q1      = imag(sigCarr .* signal1);
        I2      = real(sigCarr .* signal2);
        Q2      = imag(sigCarr .* signal2);
        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom1 = fft(I1 + 1i*Q1);
        IQfreqDom2 = fft(I2 + 1i*Q2);


        %--- Make the correlation for whole frequency band (for all freq. bins)
        for frqBinIndex = 1:numberOfFrqBins 
            % skip trailing sub-bins
            if frqBinIndex == numberOfFrqBins && binIter>1 
                continue
            end
            IQfreqDomShift1 = circshift(IQfreqDom1, frqBinIndex -1);
            IQfreqDomShift2 = circshift(IQfreqDom2, frqBinIndex -1);

            %--- Multiplication in the frequency domain (correlation in time
            % domain)
            convCodeIQ1 = IQfreqDomShift1 .* caCodeFreqDom;
            convCodeIQ2 = IQfreqDomShift2 .* caCodeFreqDom;
            %--- Perform inverse DFT and store correlation results ------------
            acqRes1 = abs(ifft(convCodeIQ1));
            acqRes2 = abs(ifft(convCodeIQ2));
            %--- Obtain current maximums 
            corrPeak1 = max(acqRes1);
            corrPeak2 = max(acqRes2);
            %--- Just keep the best correlation vector
            if corrPeak1>prevmax || corrPeak2 > prevmax
                if corrPeak1> corrPeak2
                    prevmax = corrPeak1;
                    freqShift = binIter;
                    corrVec = acqRes1;
                else
                    prevmax = corrPeak2;
                    freqShift = binIter;
                    corrVec = acqRes2;
                end
                frequencyBinIndex = frqBinIndex;
            end

        end % frqBinIndex = 1:numberOfFrqBins
    end
    % Obtain frequency bin and codePhase
    [maxPeak, codePhase] = max(corrVec);
        
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;
    
    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
            ((samplesPerBlock/Nblocks) + excludeRangeIndex1);
        
    elseif excludeRangeIndex2 >= (samplesPerBlock/Nblocks)
        codePhaseRange = (excludeRangeIndex2 - (samplesPerBlock/Nblocks) + 1) : ...
            excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
            excludeRangeIndex2 : (samplesPerBlock/Nblocks)];
    end
    
    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(corrVec(codePhaseRange));
    
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = maxPeak/secondPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (maxPeak/secondPeakSize) > settings.acqThreshold
                   
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        % Assign the obtained code phase to results
        acqResults.codePhase(PRN) = codePhase;
        % Assign the obtained Carrier Frequency to results
        acqResults.carrFreq(PRN) = initFreq - freqResolution*(frequencyBinIndex - 1) + (freqResolution/Nshifts)*(freqShift-1);
                
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
    acqResults.timeVec(idx) = toc(tstart3);
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
