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
% Modified for BDS-3 B1C by Yafeng Li
%
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $

%% Initialize constants ===================================================
numOfSatellites = size(prnList, 2);
% BDS constants
bdsPi        = 3.1415926535898;        % Pi used in BDS coordinates
%--- Constants for satellite position calculation -------------------------
OmegaE       = 7.2921150e-5;           % Earth rotation rate, [rad/s]
mu           = 3.986004418e14;         % Earth's universal gravitational
                                       % parameter,[m^3/s^2]
F            = -4.44280730904398e-10;  % Constant, [sec/(meter)^(1/2)]
A_ref_MEO    = 27906100;      % Semi-major axis reference for MEO, [m]
A_ref_IGSO_GEO = 42162200;    % Semi-major axis reference for IGSO/GEO, [m]

%% Initialize results ===============================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    
    %% Find initial satellite clock correction --------------------------------
    
    %--- Find time difference ---------------------------------------------
    dt = check_t(transmitTime(satNr) - eph(prn).t_oc);
    
    %--- Calculate clock correction ---------------------------------------
    % As joint tracking of pilot and data channel is employed, the
    % eph.ISC_B1Cd is ignored.
    satClkCorr(satNr) = (eph(prn).a_2 * dt + eph(prn).a_1) * dt + ...
        eph(prn).a_0 - ...
        eph(prn).T_GDB1Cp;
    
    time = transmitTime(satNr) - satClkCorr(satNr);
    
    %% Find satellite's position ----------------------------------------------
    % Time correction
    tk = check_t(time - eph(prn).t_oe);
    
    if strcmp(eph(prn).SatType,'MEO')
        A_ref = A_ref_MEO;
    elseif strcmp(eph(prn).SatType,'IGSO') || strcmp(eph(prn).SatType,'GEO')  
        A_ref = A_ref_IGSO_GEO;
    end

    % Restore semi-major axis
    A0   = A_ref + eph(prn).deltaA;
    A = A0 + eph(prn).ADot * tk;
    
    % Reference mean motion
    n0  = sqrt(mu / (A0^3));
    delta_n = eph(prn).delta_n_0 + 0.5 * eph(prn).delta_n_0Dot * tk;

    % Mean motion
    n = n0 + delta_n;

    % Mean anomaly
    M = eph(prn).M_0 + n * tk;
    
    % Reduce mean anomaly to between 0 and 360 deg
    M = rem(M + 2*bdsPi, 2*bdsPi);
    
    %Initial guess of eccentric anomaly
    E = M;
    
    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*bdsPi);
        
        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end
    
    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*bdsPi, 2*bdsPi);
    
    %Compute relativistic correction term
    dtr = F * eph(prn).e * sqrt(A0) * sin(E);
    
    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);
    
    %Compute angle Phi
    Phi = nu + eph(prn).omega;
    
    %Reduce phi to between 0 and 360 deg
    Phi = rem(Phi, 2*bdsPi);
    
    %Correct argument of latitude
    u = Phi + ...
        eph(prn).C_uc * cos(2*Phi) + ...
        eph(prn).C_us * sin(2*Phi);
    %Correct radius
    r = A * (1 - eph(prn).e*cos(E)) + ...
        eph(prn).C_rc * cos(2*Phi) + ...
        eph(prn).C_rs * sin(2*Phi);
    %Correct inclination
    i = eph(prn).i_0 + eph(prn).i_0Dot * tk + ...
        eph(prn).C_ic * cos(2*Phi) + ...
        eph(prn).C_is * sin(2*Phi);
    
    %Compute the angle between the ascending node and the Greenwich meridian
    Omega = eph(prn).omega_0 + (eph(prn).omegaDot - OmegaE)*tk - ...
        OmegaE * eph(prn).t_oe;
    %Reduce to between 0 and 360 deg
    Omega = rem(Omega + 2*bdsPi, 2*bdsPi);
    
    %--- Compute satellite coordinates ------------------------------------
    xp = r * cos(u);
    yp = r * sin(u);
    satPositions(1, satNr) = xp * cos(Omega) - yp * cos(i) * sin(Omega);
    satPositions(2, satNr) = xp * sin(Omega) + yp * cos(i) * cos(Omega);
    satPositions(3, satNr) = yp * sin(i);
    
    
    %% Include relativistic correction and iono in clock correction -----------
    satClkCorr(satNr) = (eph(prn).a_2 * dt + eph(prn).a_1) * dt + ...
                               eph(prn).a_0 - eph(prn).T_GDB1Cp + dtr;
    
end % for satNr = 1 : numOfSatellites
