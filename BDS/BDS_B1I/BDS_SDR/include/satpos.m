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
%                         CU Multi-GNSS SDR  
% (C) Updated by Daehee Won, Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original work by Darius Plausinaitis,Peter Rinder, 
% Nicolaj Bertelsen and Dennis M. Akos
%--------------------------------------------------------------------------


%% Initialize constants ===================================================
numOfSatellites = size(prnList, 2);

% GPS constatns

BeiDouPi       = 3.1415926535898;  % Pi used in the GPS coordinate 
                                   % system

%--- Constants for satellite position calculation -------------------------
Omegae_dot     = 7.2921150e-5;     % Earth rotation rate, [rad/s]
GM             = 3.986004418e14;   % Earth's universal gravitational constant
                                   % gravitational parameter,
                                   % [m^3/s^2]
F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    
%% Find initial satellite clock correction --------------------------------

    %--- Find time difference ---------------------------------------------
%     dt = check_t(transmitTime(satNr) - eph(prn).t_oc);
    dt = check_t(transmitTime(satNr) - eph(prn).t_oc);

    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph(prn).a2 * dt + eph(prn).a1) * dt + ...
                         eph(prn).a0 - ...
                         eph(prn).T_GD_1;

%     time = transmitTime(satNr) - satClkCorr(satNr);
    time = transmitTime(satNr) - satClkCorr(satNr);

%% Find satellite's position ----------------------------------------------

    %Restore semi-major axis
    A   = eph(prn).sqrtA * eph(prn).sqrtA;

    %Time correction
    tk  = check_t(time - eph(prn).t_oe);

    %Initial mean motion
    n0  = sqrt(GM / A^3);
    %Mean motion
    n   = n0 + eph(prn).deltan;

    %Mean anomaly
    M   = eph(prn).M_0 + n * tk;
    %Reduce mean anomaly to between 0 and 360 deg
    M   = rem(M + 2*BeiDouPi, 2*BeiDouPi);

    %Initial guess of eccentric anomaly
    E   = M;

    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*BeiDouPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end

    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*BeiDouPi, 2*BeiDouPi);

    %Compute relativistic correction term
    dtr = F * eph(prn).e * eph(prn).sqrtA * sin(E);
%     tk  = check_t(time + dtr - eph(prn).t_oe);    

    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);

    %Compute angle phi
    phi = nu + eph(prn).omega;
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*BeiDouPi);

    %Correct argument of latitude
    u = phi + ...
        eph(prn).C_uc * cos(2*phi) + ...
        eph(prn).C_us * sin(2*phi);
    %Correct radius
    r = A * (1 - eph(prn).e*cos(E)) + ...
        eph(prn).C_rc * cos(2*phi) + ...
        eph(prn).C_rs * sin(2*phi);
    %Correct inclination
    i = eph(prn).i_0 + eph(prn).iDot * tk + ...
        eph(prn).C_ic * cos(2*phi) + ...
        eph(prn).C_is * sin(2*phi);

    
    %% GEO Satellite coordinates
    if prn <= 5
                
        %Compute the angle between the ascending node and the Greenwich meridian
        Omega = eph(prn).omega_0 + eph(prn).omegaDot*tk - ...
                Omegae_dot * eph(prn).t_oe;
        %Reduce to between 0 and 360 deg
        Omega = rem(Omega + 2*BeiDouPi, 2*BeiDouPi);
        
        %--- Compute satellite coordinates ------------------------------------
        satPositions(1, satNr) = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
        satPositions(2, satNr) = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
        satPositions(3, satNr) = sin(u)*r * sin(i);
        
        angX = -5*(BeiDouPi/180);
        angZ = Omegae_dot*tk;
        
        Rx = [1      0         0;
              0  cos(angX) sin(angX);
              0 -sin(angX) cos(angX)];
        
        Rz = [cos(angZ) sin(angZ) 0;
             -sin(angZ) cos(angZ) 0;
                  0         0     1];
        satPositions(:, satNr) = Rz * Rx * satPositions(:, satNr);
        
        
    %% MEO/IGSO Satellite coordinates    
    elseif prn > 5
        
        %Compute the angle between the ascending node and the Greenwich meridian
        Omega = eph(prn).omega_0 + (eph(prn).omegaDot - Omegae_dot)*tk - ...
                Omegae_dot * eph(prn).t_oe;
        %Reduce to between 0 and 360 deg
        Omega = rem(Omega + 2*BeiDouPi, 2*BeiDouPi);
        
        %--- Compute satellite coordinates ------------------------------------
        satPositions(1, satNr) = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
        satPositions(2, satNr) = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
        satPositions(3, satNr) = sin(u)*r * sin(i);
    end


%% Include relativistic correction in clock correction --------------------
    satClkCorr(satNr) = (eph(prn).a2 * dt + eph(prn).a1) * dt + ...
                         eph(prn).a0 - ...
                         eph(prn).T_GD_1 + dtr;
                     
end % for satNr = 1 : numOfSatellites
