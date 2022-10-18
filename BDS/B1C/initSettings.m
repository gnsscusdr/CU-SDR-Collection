function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".  
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
% and Dennis M. Akos. 
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos.
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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Raw signal file name and related parameter =======================
settings.fileName           = '/media/gnss/Ext2TB/data/L1/B200Ex_1.bin';
% Data type used to store one sample
settings.dataType           = 'schar';
% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...                      
settings.fileType           = 2;
%Intermediate, sampling, code freqency and code length 
settings.IF                 = -20e3;  %[Hz]
settings.samplingFreq       = 18e6;              %[Hz]
% Front-end bandwidth for calculation of weighting factor in wideband (WB)
% tracking mode (this is required only for WB tracking mode)
settings.FEBW   = 27e6;          % [Hz]

%% Important operation modes configuration ==========================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 40000;        %[ms]
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition. In May 2020, BeiDou3 includes PRNs with 19:30,32:46,59,60.
settings.acqSatelliteList   = [1:69]; 
% Enable use of pilot channel to do acquisition
settings.pilotACQflag       = 1;             % 0 - Off; 1 - On
% Pilot channel tracking mode setting
settings.pilotTRKflag       = 1;             % 1 - Narrowband; 2 - Full band 
% Number of channels to be used for signal processing
settings.numberOfChannels   = 15;
% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
settings.skipNumberOfBytes  = 0e6;
                                             
%% B1C signal parameters ============================================
% Define number of chips in a B1C primary code period.
settings.codeLength    = 10230;        % [Chips]
% Nominal rate of B1C primary code. 
settings.codeFreqBasis = 1.023e6;      % [Hz]
% Carrier center frequency
settings.carrFreqBasis = 1575.42e6;    % [Hz]

%% Acquisition settings =============================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% It is single sideband, so the whole search band is tiwce of it.
settings.acqSearchBand      = 5000;           % [Hz]
% Coherent integration time for acqisition: sub-portion of code used to do
% correlation. This parameter should be greater than 1ms and less than 10 ms
settings.acqCohT            = 10;             % [ms]
% Frequency search step
settings.acqStep            = 1000/settings.acqCohT/2;  % [Hz] 
% Threshold for the signal presence decision rule
settings.acqThreshold       = 10;
% Sampling rate threshold for downsampling 
settings.resamplingThreshold    = 15e6;           % [Hz]
% Enable/dissable use of downsampling for acquisition
settings.resamplingflag         = 0;              % 0 - Off; 1 - On    
                                                  
%% Tracking loops settings ==========================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 1;      % [Hz]
% Hailf of the E-L spacing should be kept below approximately 0.0765 chip
% for full-band tracking mode (see Betz 2016)
settings.dllCorrelatorSpacing    = 0.06;     % [primary code chip] 
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 18;       % [Hz]
% Precorrelation integration time for DLL and PLL
settings.intTime = 0.01;                     % [s]

%% Navigation solution settings =====================================
% Period for calculating pseudoranges and position
settings.navSolPeriod       = 200;          %[ms]
% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 5;           %[degrees 0 - 90]
% Enable/disable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off; 1 - On
% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking    = 1;            % 0 - Off; 1 - On
settings.plotAcquisition = 1;
settings.plotNavigation  = 1;

%% Constants ========================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       % [ms] Initial sign. travel time

%% CNo Settings======================================================
% Number of correlation values used to compute each C/No point
settings.CNoInterval = 50;
