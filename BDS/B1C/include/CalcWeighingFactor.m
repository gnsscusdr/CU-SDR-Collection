function factor = CalcWeighingFactor(settings)
%This function calculates the weighting factor of data channel for 
%constructing the composite code tracking error in fullband tracking mode.  
%
%factor = weighingFactor(settings)
%
%   Inputs:
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       factor       - weighting factor of data channel 
 
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos

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

%CVS record:
%$Id: calcLoopCoef.m,v 1.1.2.2 2018/06/25 11:38:22 dpl Exp $

% Code frequency and code chip length -------------------------------------
fc = settings.codeFreqBasis;
Tc = 1/fc;
% Front end bandwidth
Br = settings.FEBW;

% Power ratio and RMS BW of data channel ----------------------------------
% PSD of data channel: BOC(1,1)
G_BOC1_1f = @(f) Tc .* ( sin(pi/2*f/fc).* sin(pi*f/fc)./cos(pi/2*f/fc)*fc./f/pi ).^2;
% For calculating RMS BW of data channel
G_BOC1_1f_2 = @(f) (Tc .* ( sin(pi/2*f/fc).* sin(pi*f/fc)./cos(pi/2*f/fc)*fc./f/pi ).^2 .*f.^2);

% For calculating RMS BW of data channel
Power_BOC1_1_2 = integral(@(f)G_BOC1_1f_2(f),-Br/2,Br/2);

% Power ratio of data channel 
Power_BOC1_1 = integral(@(f)G_BOC1_1f(f),-Br/2,Br/2);
% RMS BW of data channel
REM_BW_BOC_1_1 = (Power_BOC1_1_2 ./ Power_BOC1_1).^0.5;

% Power ratio and RMS BW of pilot channel ---------------------------------
% PSD of pilot channel: QMBOC
G_pilotf = @(f) 29/33*Tc .* ( sin(pi/2*f/fc).* sin(pi*f/fc)./cos(pi/2*f/fc)*fc./f/pi ).^2 + ...
    4/33*Tc .* ( sin(pi/12*f/fc).* sin(pi*f/fc)./cos(pi/12*f/fc)*fc./f/pi ).^2;
% For calculating RMS BW of pilot channel
G_pilotf_2 = @(f) (29/33*Tc .* ( sin(pi/2*f/fc).* sin(pi*f/fc)./cos(pi/2*f/fc)*fc./f/pi ).^2 .*f.^2 + ...
    4/33*Tc .* ( sin(pi/12*f/fc).* sin(pi*f/fc)./cos(pi/12*f/fc)*fc./f/pi ).^2 .*f.^2) ;

% For calculating RMS BW of pilot channel
Power_pilot_2 = integral(@(f)G_pilotf_2(f),-Br/2,Br/2);

% Power ratio of pilot channel
Power_pilot = integral(@(f)G_pilotf(f),-Br/2,Br/2);
% RMS BW of pilot channel
REM_BW_pilot = (Power_pilot_2 ./ Power_pilot).^0.5;

% Calculate the weighting factor ------------------------------------------
temp1 = 11 * Power_BOC1_1 * REM_BW_BOC_1_1^2;
temp2  = 33 * Power_pilot * REM_BW_pilot^2;
factor =  temp1/ (temp1 + temp2);

