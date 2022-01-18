function B3ICodesTable = makeB3ITable(PRN,settings)
%Function generates B3I codes for all specified PRN based on the settings
%provided in the structure "settings". The codes are digitized at the
%sampling frequency specified in the settings structure.
%
%B3ICodesTable = makeB3ITable(PRN,settings)
%
%   Inputs:
%       PRN             - specified PRN for B3I code
%       settings        - receiver settings
%   Outputs:
%       B3ICodesTable   - an vector containing the sampled B3I codes
%                       for specified PRN

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B3I SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
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
%$Id: makeCaTable.m,v 1.1.2.6 2020/01/16 11:38:22 dpl Exp $

%--------------------------------------------------------------------------
% Modified for Beidou by Yafeng Li
% Final update: 2020/01/20

%--- Find number of samples per spreading code ----------------------------
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

%--- Find time constants --------------------------------------------------
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % B3I chip period in sec

%--- Generate B3I code for given PRN --------------------------------------
B3ICode = generateB3Icode(PRN);

%=== Digitizing ===========================================================
%--- Make index vector to read B3I code values ----------------------------
% The length of the index vector depends on the sampling frequency -
% number of samples per millisecond (because one B3I code period is one
% millisecond).
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

%--- Correct the last index (due to number rounding issues) ---------------
codeValueIndex(end) = 10230;  % B3I chip lenth is 10230

%--- Make the digitized version of the B3I code ---------------------------
% The "upsampled" code is made by selecting values form the B3I code
% chip vector (B3ICode) for the time instances of each sample.
B3ICodesTable = B3ICode(codeValueIndex);
