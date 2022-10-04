% Script postProcessing.m processes the raw signal from the specified data
% file (in settings) operating on blocks of 37 seconds of data.
%
% First it runs acquisition code identifying the satellites in the file,
% then the code and carrier for each of the satellites are tracked, storing
% the 1msec accumulations.  After processing all satellites in the 37 sec
% data block, then postNavigation is called. It calculates pseudoranges
% and attempts a position solutions. At the end plots are made for that
% block of data.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Updated by Jakob Almqvist, Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original work by Darius Plausinaitis,Peter Rinder, 
% Nicolaj Bertelsen and Dennis M. Akos
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

%                         THE SCRIPT "RECIPE"
%
% The purpose of this script is to combine all parts of the software
% receiver.
%
% 1.1) Open the data file for the processing and seek to desired point.
%
% 2.1) Acquire satellites
%
% 3.1) Initialize channels (preRun.m).
% 3.2) Pass the channel structure and the file identifier to the tracking
% function. It will read and process the data. The tracking results are
% stored in the trkResults structure. The results can be accessed this
% way (the results are stored each millisecond):
% trkResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
% XXX is a field name of the result (e.g. I_P, codePhase etc.)
%
% 4) Pass tracking results to the navigation solution function. It will
% decode navigation messages, find satellite positions, measure
% pseudoranges and find receiver position.
%
% 5) Plot the results.

%% Initialization =========================================================
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
    
    % Find number of samples per spreading code
    samplesPerCode = round(settings.samplingFreq / ...
                       (settings.codeFreqBasis / settings.codeLength));
    
    
    % At least 42ms of signal are needed for fine frequency estimation
    codeLen = max(42,settings.acqNonCohTime+2);
    % Read data for acquisition.
    data  = fread(fid, dataAdaptCoeff*codeLen*samplesPerCode, settings.dataType)';

    if (dataAdaptCoeff==2)    
        data1=data(1:2:end);    
        data2=data(2:2:end);    
        data=data1 + 1i .* data2;    
    end

    %--- Do the acquisition -------------------------------------------
    disp ('   Acquiring satellites...');
    acqResults = acquisition(data, settings);
    save("acqResults")
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
[trkResults, ~] = tracking(fid, channel, settings);
save("trkResults")
% Close the data file
fclose(fid);
disp(['   Tracking is over (elapsed time ', ...
                                    datestr(now - startTime, 13), ')'])                      

%% Calculate navigation solutions =========================================
disp('   Calculating navigation solutions...');

[navResults, ~] = postNavigation(trkResults, settings);
save("navResults")

disp('   Processing is complete for this data block');
disp('Post processing of the signal is over.');
%% Plot all results ===================================================
disp ('   Ploting results...');

if settings.plotAcquisition
    plotAcquisition(acqResults);
end

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

