function PilotCodesTable = makePilotTable(settings,PRN)
%Function generates BDS-3 B1C BOC(1,1) codes of data chnnel for specified
%PRN based on the settings provided in the structure "settings". The codes
%are digitized at the sampling frequency specified in the settings structure.
%
%PilotCodesTable = makePilotTable(settings,PRN)
%
%   Inputs:
%       settings         - receiver settings
%   Outputs:
%       PilotCodesTable  - sampled BDS-3 B1C BOC(1,1) spreading waveform
%                          for pilot channel

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Developed for BDS B1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
% and Dennis M. Akos. 
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos

% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and 
% implementation of an open-source BDS-3 B1C/B2a SDR receiver. 
% GPS Solut 23, 60 (2019). https://doi.org/10.1007/s10291-019-0853-z
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
%$Id: makeCaTable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;       % Sampling period in sec
tc = 1/settings.codeFreqBasis/2;    % B1C BOC(1,1) chip period in sec


%--- Generate B1C code for given PRN -----------------------------------
PilotCode = generatePilotBOC11(settings,PRN);

%=== Digitizing =======================================================

%--- Make index array to read B1C code values -------------------------
% The length of the index array depends on the sampling frequency
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

%--- Correct the last index (due to number rounding issues) -----------
codeValueIndex(end) = settings.codeLength*2;
codeValueIndex(1) = 1;

%--- Make the digitized version of the B1C code -----------------------
% The "upsampled" code is made by selecting values from the B1C code
% chip array for the time instances of each sample.
PilotCodesTable = PilotCode(codeValueIndex);



