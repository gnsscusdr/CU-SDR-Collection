function [pf3,pf2,pf1] = calcLoopCoefCarr(settings)
%Function finds loop coefficients for 3-order loop. The coefficients are 
% used then in PLL-s to improve the dynamic performance of carrier tracking.
%
%[pf3,pf2,pf1] = calcLoopCoefCarr(settings)
%
%   Inputs:
%       settings      - receiver settings.
%
%   Outputs:
%       pf3,pf2,pf1   - Three Loop filter coefficients 
 
%--------------------------------------------------------------------------
% Written by Yafeng Li. 
% See reference "Satellite Signal Acquisition, Tracking, and Data 
% Demodulation" by Kaplan and Hegarty.
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

% Loop noise bandwidth: The 3-order loop remains stable at LBW <= 18Hz.  
LBW = settings.pllNoiseBandwidth;

% Summation interval
intTime = settings.intTime;

% loop constant coefficients
a3 = 1.1;
b3 = 2.4;

% Solve natural frequency
Wn = LBW/0.7845;

% solve for [pf3,pf2,pf1]
pf3 = Wn^3 * intTime^2;
pf2 = a3 * Wn^2 * intTime;
pf1 = b3 * Wn;
