disp ('Starting processing...');

[fid, message] = fopen(settings.fileName, 'rb');

%Initialize the multiplier to adjust for the data type
if (settings.fileType==1)
dataAdaptCoeff=1;
else
dataAdaptCoeff=2;
end

%If success, then process the data
if (fid > 0)

% Move the starting point of processing. Can be used to start the
% signal processing at any point in the data record (e.g. good for long
% records or for signal processing in blocks).
fseek(fid, dataAdaptCoeff*settings.skipNumberOfBytes, 'bof');

%% Acquisition ============================================================

% Do acquisition if it is not disabled in settings or if the variable
% acqResults does not exist.
if ((settings.skipAcquisition == 0) || ~exist('acqResults', 'var'))
    %--- Read data for acquisition ------------------------------------
    % Find number of samples per spreading code
    samplesPerCode = round(settings.samplingFreq / ...
        (settings.codeFreqBasis / settings.codeLength));
    
    % Read data for acquisition: use 1 code length data to do fine
    % acquisition; To ensure there is at least 1 complete B1C code is
    % included, 2 code length of data is needed
    data  = fread(fid, dataAdaptCoeff*2*samplesPerCode, settings.dataType)';
    
    % Covert complex data into normal representation.
    if (dataAdaptCoeff==2)
        data1=data(1:2:end);
        data2=data(2:2:end);
        data=data1 + 1i .* data2;
    end
    
    %--- Do the acquisition -------------------------------------------
    disp ('   Acquiring satellites...');
    % Normal acquisition
    acqResults = acquisitionCS(data, settings);        
end

%% Initialize channels and prepare for the run ============================

% Start further processing only if a GNSS signal was acquired (the
% field FREQUENCY will be set to 0 for all not acquired signals)
if (any(acqResults.carrFreq))
    channel = preRun(acqResults, settings);
    showChannelStatus(channel, settings);
else
    % No satellites to track, exit
    disp('No GNSS signals detected, signal processing finished.');
    trkResults = [];
    navResults = [];
    return;
end

%% Track the signal =======================================================
startTime = now;
disp (['   Tracking started at ', datestr(startTime)]);

% Process all channels for given data block
if settings.pilotTRKflag == 1     % Narrowband tracking mode
    % Tracking with Matlab correlator
    [trkResults, ~] = NB_tracking(fid, channel, settings);
elseif settings.pilotTRKflag == 2  % Wideband tracking mode
    % Tracking with Matlab correlator
    [trkResults, ~] = WB_tracking(fid, channel, settings);
end
% Close the data file
fclose(fid);

disp(['   Tracking is over (elapsed time ', ...
    datestr(now - startTime, 13), ')'])


%% Calculate navigation solutions =========================================
disp('   Calculating navigation solutions...');
[navResults, ~] = postNavigation(trkResults, settings);
disp('   Processing is complete for this data block');

%% Plot all results ===================================================
disp ('   Ploting results...');
if settings.plotTracking
    plotTracking(1:settings.numberOfChannels, trkResults, settings);
end

if settings.plotNavigation
    plotNavigation(navResults, settings);
end
disp('Post processing of the signal is over.');

else
% Error while opening the data file.
error('Unable to read file %s: %s.', settings.fileName, message);
end % if (fid > 0)


