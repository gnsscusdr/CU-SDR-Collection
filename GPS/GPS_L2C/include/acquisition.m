function acqResults = acquisition(longSignal, settings)
%% Acquisition initialization ===============================================
%--- Find number of samples for fiffernet long signals ------------------------------
Nblocks = 2;
% Number of samples per L2CM code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Number of samples per chip
samplesPerChip = round(settings.samplingFreq / settings.codeFreqBasis);
% Samples per block
samplesPerBlock = samplesPerCode*Nblocks;
%--- Cut 20 plus X cm input signal to do correlate ----------------------------------
signal = longSignal(1:samplesPerBlock);

%--- Generate input and local signal to to correlate ------------------------
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of the local carrier wave
phasePoints = (0 : (samplesPerBlock-1)) * 2 * pi * ts;
% Obtain frequency resolution
freqResolution = settings.samplingFreq/samplesPerBlock;
% Number of the frequency bins for the given acquisition band
numberOfFrqBins = round(settings.acqSearchBand*1e3/freqResolution) + 1;
% Determine number of sub-bins
Nshifts = freqResolution/settings.acqStep;
%% Initialize variables
%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);
% Initialize the starting frequency
initFreq = settings.IF + (settings.acqSearchBand/2) * 1000;

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList    
    % Initialize variables
    prevmax = 0;
    corrVec = zeros(1, samplesPerBlock);
    frequencyBinIndex = 0;
    results = zeros(numberOfFrqBins*Nshifts - Nshifts, samplesPerBlock);
    freqShift = 0;
    % Generate CM code and local replica for current PRN
    cmCodesTable = makeCMTable(settings,PRN);
    localCode = [cmCodesTable(1:samplesPerCode), zeros(1,samplesPerCode)];
    cmCodeFreqDom = conj(fft(localCode));

    for binIter = 1:Nshifts
        %--- Generate carrier wave frequency grid (1kHz step per 1msec) -----------
        initFreqShift = initFreq - (binIter-1)*(freqResolution/Nshifts);
        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(-1i*initFreqShift*phasePoints);
        %--- "Remove carrier" from the signal -----------------------------
        I1      = real(sigCarr .* signal);
        Q1      = imag(sigCarr .* signal);
        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom = fft(I1 + 1i*Q1);
    
        %% Circular shift
        %--- Make the correlation for whole frequency band (for all freq. bins)
        for frqBinIndex = 1:numberOfFrqBins
            % skip trailing sub-bins to avoid problems with matrix dimensions
            if frqBinIndex == numberOfFrqBins && binIter>1 
                continue
            end
            % Shift local replica
            IQfreqDomShift = circshift(IQfreqDom, frqBinIndex -1);
            %--- Multiplication in the frequency domain (correlation in time
            %domain)
            convCodeIQ = IQfreqDomShift .* cmCodeFreqDom;
            %--- Perform inverse DFT and store correlation results ------------
            acqRes = abs(ifft(convCodeIQ)); 
            % Obtain current maximum
            currmax = max(acqRes);
            if currmax > prevmax
                prevmax = currmax;
                corrVec = acqRes;
                frequencyBinIndex = frqBinIndex;
                freqShift = binIter;
            end  
        end % frqBinIndex = 1:numberOfFrqBins

    end
    % Obtain codePhase and Peak
    [maxPeak,codePhase] = max(corrVec);
    %--- Find 1 chip wide CM code phase exclude range around the peak ----
    excludeRangeIndex1 = codePhase - samplesPerChip;
    excludeRangeIndex2 = codePhase + samplesPerChip;
    
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
        acqResults.carrFreq(PRN) = initFreq - freqResolution*(frequencyBinIndex - 1) - (freqResolution/Nshifts)*(freqShift-1);
        acqResults.codePhase(PRN) = codePhase;

        % ========== Find the L2CL code phase ===========         
        % Using detected CM code phase
        signal0DC = longSignal(codePhase:(codePhase + samplesPerCode-1));
        signal0DC = signal0DC - mean(signal0DC);

        if (settings.pilotTRKflag == 1)
            
            powerArray = zeros(1,75);
            
            % Find phase points of the local carrier wave
            phasePointsCL = (0 : (samplesPerCode-1)) * 2 * pi * ts;
            
            % Generate local sine and cosine
            sigCarr = exp(-1i* acqResults.carrFreq(PRN) * phasePointsCL);
            
            CLCode = generateCLcode(PRN,settings);
            
            %--- Find time constants --------------------------------------------------
            tc = 1/(settings.codeFreqBasis*2);   % PRN chip period in sec
            
            %=== Digitizing =======================================================
            
            %--- Make index array to read C/A code values -------------------------
            % The length of the index array depends on the sampling frequency -
            % number of samples per millisecond (because one C/A code period is one
            % millisecond).
            codeValueIndex = ceil((ts * (0:samplesPerCode-1)) / tc);
            
            %--- Correct the last index (due to number rounding issues) -----------
            codeValueIndex(1) = 1;
            if(settings.acqCohT <= 10)
                % it is settings.codeLength/20 * 10 * 2
                codeValueIndex(end) = settings.codeLength;
            else
                codeValueIndex(end) = settings.codeLength *2;
            end
            
            for ind = 1:75
                %--- Make the digitized version of the CL code -----------------------
                % The "upsampled" code is made by selecting values form the CA code
                % chip array (caCode) for the time instances of each sample.
                CLCodeSample = CLCode(codeValueIndex + settings.codeLength * 2* (ind-1));
                
                powerArray(1,ind) = abs(sum(signal0DC.*CLCodeSample.*sigCarr));
            end
            [~,acqResults.CLCodePhase(PRN)]  = max(powerArray);
        end
            
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
