function [CNo, PllDetector]= Calc_CNo_PLD(trackResults,settings,loopCnt)
%Calculate CNo using the Variance Summing Method, and the PLL lock detector
%output
%
%[CNo, PllDetector]= Calc_CNo_PLD(trackResults,settings,loopCnt)
%
%   Inputs:
%       trackResults      - Correlation values
%       settings          - Receiver settings
%       loopCnt           - Iteration index for C/No calculation
%   Outputs:
%       CNo               - Estimated C/No for the given values of I and Q
%                           CNo(1) is for data channel,CNo(21) is for pilot
%                           channel, and CNo(2) is for whole B1C signal  
%       PllDetector       - PLL lock detector output for data and pilot
%                            channels
%
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

CNo = zeros(1,3);
PllDetector = zeros(1,2);
% Integration time for each correlation point
T = settings.intTime;
%% ---- C/No and PLL detector output estimations for data channel -----
% Data-channel prompt correlation values
I_P = trackResults.I_P(loopCnt - settings.CNoInterval + 1 : loopCnt);
Q_P = trackResults.Q_P(loopCnt - settings.CNoInterval + 1 : loopCnt);

% ----------------------- CNo extimation for data channel -----------------
% Calculate Power
Z = I_P.^2 + Q_P.^2;
% Calculate the mean and variance of the Power
Zm = mean(Z);
Zv = var(Z);
% Calculate the average carrier power
Pav = sqrt(Zm^2 - Zv);
% Calculate the variance of the noise
Nv = 0.5 * (Zm - Pav);
%  C/No estimation for data channel
DataCNo = abs((1/T) * Pav / (2 * Nv));
CNo(1) = 10*log10(DataCNo);

% ----------------------- PLL lock detector output -----------------------
% Narrowbnd power
NBP = (sum(I_P(I_P>0)) - sum(I_P(I_P<0)))^2 + sum(Q_P)^2;
NBD = (sum(I_P(I_P>0)) - sum(I_P(I_P<0)))^2 - sum(Q_P)^2;
% PLL lock detector output
PllDetector(1) = NBD/NBP;

%% ---- C/No and PLL detector output estimations for pilot channel ----

% CNo for pilot channel in unit of ratio-Hz
PilotCNo = 0;

if settings.pilotTRKflag == 2
    % Pilot-channel prompt correlation values
    I_P = trackResults.Pilot_I_P(loopCnt - settings.CNoInterval + 1 : loopCnt);
    Q_P = trackResults.Pilot_Q_P(loopCnt - settings.CNoInterval + 1 : loopCnt);
elseif settings.pilotTRKflag == 1
    % Pilot-channel prompt correlation values
    Q_P = trackResults.Pilot_I_P(loopCnt - settings.CNoInterval + 1 : loopCnt);
    I_P = trackResults.Pilot_Q_P(loopCnt - settings.CNoInterval + 1 : loopCnt);
end

if (settings.pilotTRKflag == 1) || (settings.pilotTRKflag == 2)    
    % -------------------- CNo extimation for pilot channel ---------------
    % Calculate Power
    Z = I_P.^2 + Q_P.^2;
    % Calculate the mean and variance of the Power
    Zm = mean(Z);
    Zv = var(Z);
    % Calculate the average carrier power
    Pav = sqrt(Zm^2 - Zv);
    % Calculate the variance of the noise
    Nv = 0.5 * (Zm - Pav);
    %  C/No estimation for data channel
    PilotCNo = abs((1/T) * Pav / (2 * Nv));
    CNo(2) = 10*log10(PilotCNo);
    
    % ------------ PLL lock detector output for pilot channel -------------
    % Narrowbnd power
    NBP = (sum(I_P(I_P>0)) - sum(I_P(I_P<0)))^2 + sum(Q_P)^2;
    NBD = (sum(I_P(I_P>0)) - sum(I_P(I_P<0)))^2 - sum(Q_P)^2;
    % PLL lock detector output
    PllDetector(2) = NBD/NBP;
end

% CNo extimation for total B1C siganl channel
CNo(3) = 10 * log10(DataCNo + PilotCNo );
