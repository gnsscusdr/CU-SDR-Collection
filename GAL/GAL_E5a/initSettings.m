function settings = initSettings()

%% Processing settings ====================================================
% Number of milliseconds to be processed used 32000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 65000;        %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels   = 12;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
settings.skipNumberOfBytes     = 0;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode

settings.fileName           = '../../../dataSets/L5_IF20KHz_FS18MHz.bin';

% Data type used to store one sample
settings.dataType           = 'schar';

% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...                      
settings.fileType           = 2;

% Intermediate, sampling and code frequencies
settings.IF                 = -20e3;        % [Hz] (Invert sign here for proper acquisition)
settings.samplingFreq       = 18e6;         % [Hz]  
settings.codeFreqBasis      = 10.23e6;      % [Hz]
% Define number of chips in a E5aI/E5aQ tiered code period
settings.codeLength          = 10230;

%% Acquisition settings ===================================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition. As of June 2020, in-orbit Galileo SVs includes PRNs:
% [1 2 3 4 5 7 8 9 11 12 13 14 15 18 19 21 22 24 25 26 27 30 31 33 36]
settings.acqSatelliteList   = [1:36]; 
% Band around IF to search for satellite signal. Depends on max Doppler.
% It is single sideband, so the whole search band is tiwce of it.
settings.acqSearchBand          = 5000;           % [Hz]
% Non-coherent integration times after 1ms coherent integration
settings.acqNonCohTime      = 15;              %[ms]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 4.5;
% Frequency search step for coarse acquisition
settings.acqSearchStep      = 500;                % [Hz]
% Sampling rate threshold for downsampling 
settings.resamplingThreshold    = 45e6;           % [Hz]
% Enable/dissable use of downsampling for acquisition
settings.resamplingflag         = 0;              % 0 - Off; 1 - On
%% Tracking loops settings ================================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 1.5;%2     %[Hz]
settings.dllCorrelatorSpacing    = 0.5;     %[chips]

% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 15;  %18     %[Hz]
% Integration time for DLL and PLL
settings.intTime = 0.001;%0.02;
% Enable/dissable use of pilot channel for tracking
settings.pilotTRKflag            = 1;        % 0 - Off
                                             % 1 - On
%% Navigation solution settings ===========================================
% Period for calculating pseudoranges and position
settings.navSolPeriod       = 500;          %[ms]

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 5;           %[degrees 0 - 90]
% Enable/disable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off
                                            % 1 - On

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;
%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off
                                            % 1 - On
settings.plotAcquisition       = 1;      
settings.plotNavigation       = 1;      
%% Constants ==============================================================

settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time

%% CNo Settings============================================================
% Accumulation interval in Tracking (in Sec)
settings.CNo.accTime = 0.001;%0.02;
% Accumulation interval for computing VSM C/No (in ms)
settings.CNo.VSMinterval = 100;%20;
%% E5b carrier frequency ===================================================
settings.carrFreqBasis    = 1176.45e6;       %[Hz]