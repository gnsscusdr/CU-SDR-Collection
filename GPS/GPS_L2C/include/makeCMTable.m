function cmCodesTable = makeCMTable(settings,PRN)
%Function generates CA codes for all 32 satellites based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%One row in the "caCodesTable" is one C/A code. The row number is the PRN
%number of the C/A code.
%
%caCodesTable = makeCaTable(settings)
%
%   Inputs:
%       settings        - receiver settings
%   Outputs:
%       caCodesTable    - an array of arrays (matrix) containing C/A codes
%                       for all satellite PRN-s

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

%CVS record:
%$Id: makeCaTable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
 
%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;        % Sampling period in sec
tc = 1/(settings.codeFreqBasis*2);   % PRN chip period in sec
 
%--- Generate CM code for given PRN -----------------------------------------
cmCode = generateCMcode(PRN,settings);

%=== Digitizing =======================================================

%--- Make index array to read C/A code values -------------------------
% The length of the index array depends on the sampling frequency -
% number of samples per millisecond (because one C/A code period is one
% millisecond).
codeValueIndex = ceil((ts * (0:samplesPerCode-1)) / tc);

%--- Correct the last index (due to number rounding issues) -----------
codeValueIndex(end) = settings.codeLength*2;
codeValueIndex(1) = 1;

%--- Make the digitized version of the C/A code -----------------------
% The "upsampled" code is made by selecting values form the CA code
% chip array (caCode) for the time instances of each sample.
cmCodesTable = cmCode(codeValueIndex);
