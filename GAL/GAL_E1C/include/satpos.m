function [satPositions, satClkCorr] = satpos(transmitTime, prnList, ...
                                             eph) 
%SATPOS Calculation of X,Y,Z satellites coordinates at TRANSMITTIME for
%given ephemeris EPH. Coordinates are calculated for each satellite in the
%list PRNLIST.
%[satPositions, satClkCorr] = satpos(transmitTime, prnList, eph, settings);
%
%   Inputs:
%       transmitTime  - transmission time
%       prnList       - list of PRN-s to be processed
%       eph           - ephemeridies of satellites
%       settings      - receiver settings
%
%   Outputs:
%       satPositions  - positions of satellites (in ECEF system [X; Y; Z;])
%       satClkCorr    - correction of satellites clocks

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre 04-09-96
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $

%% Initialize constants ===================================================
numOfSatellites = size(prnList, 2);

% Galileo constants
galPi          = 3.1415926535898;  % Pi used in Galileo coordinates 
                                   
%--- Constants for satellite position calculation -------------------------
OmegaE         = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
mu             = 3.986004418e14;      % Earth's universal
                                   % gravitational parameter,
                                   % [m^3/s^2]
F              = -4.442807309e-10; % Constant, [sec/(meter)^(1/2)]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    
%% Find initial satellite clock correction --------------------------------

    %--- Find time difference ---------------------------------------------
    dt = check_t(transmitTime(satNr) - eph(prn).t_oc);

    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - ...
                         eph(prn).BGD_E1E5b;

    time = transmitTime(satNr) - satClkCorr(satNr);

%% Find satellite's position ----------------------------------------------

    % Restore semi-major axis
    A   = (eph(prn).sqrtA)^2;
    
    % Reference mean motion
    n0  = sqrt(mu / (A^3));
    
    % Mean motion
    n = n0 + eph(prn).deltan;

    % Time correction
    tk = check_t(time - eph(prn).t_oe);
    
    % Mean anomaly
    M = eph(prn).M_0 + n * tk;
    
    % Reduce mean anomaly to between 0 and 360 deg
    M = rem(M + 2*galPi, 2*galPi);

    %Initial guess of eccentric anomaly
    E = M;

    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*galPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end

    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*galPi, 2*galPi);

    %Compute relativistic correction term
    dtr = F * eph(prn).e * eph(prn).sqrtA * sin(E);

    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);

    %Compute angle Phi
    Phi = nu + eph(prn).omega;
    
    %Reduce phi to between 0 and 360 deg
    Phi = rem(Phi, 2*galPi);

    %Correct argument of latitude
    u = Phi + ...
        eph(prn).CUC * cos(2*Phi) + ...
        eph(prn).CUS * sin(2*Phi);
    %Correct radius
    r = A * (1 - eph(prn).e*cos(E)) + ...
        eph(prn).CRC * cos(2*Phi) + ...
        eph(prn).CRS * sin(2*Phi);
    %Correct inclination
    i = eph(prn).i_0 + eph(prn).iDot * tk + ...
        eph(prn).CIC * cos(2*Phi) + ...
        eph(prn).CIS * sin(2*Phi);
    
    %Compute the angle between the ascending node and the Greenwich meridian
    Omega = eph(prn).Omega_0 + (eph(prn).OmegaDot - OmegaE)*tk - ...
            OmegaE * eph(prn).t_oe;
    %Reduce to between 0 and 360 deg
    Omega = rem(Omega + 2*galPi, 2*galPi);

    %--- Compute satellite coordinates ------------------------------------
    xp = r * cos(u);
    yp = r * sin(u);
    satPositions(1, satNr) = xp * cos(Omega) - yp * cos(i) * sin(Omega);
    satPositions(2, satNr) = xp * sin(Omega) + yp * cos(i) * cos(Omega);
    satPositions(3, satNr) = yp * sin(i);


%% Include relativistic correction and iono in clock correction -----------
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - ...
                         eph(prn).BGD_E1E5b + dtr;
                     
end % for satNr = 1 : numOfSatellites
