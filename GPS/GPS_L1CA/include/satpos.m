function [satPositions, satClkCorr] = satpos(transmitTime, prnList,eph) 
%SATPOS Calculation of X,Y,Z satellites coordinates at TRANSMITTIME for
%given ephemeris EPH. Coordinates are calculated for each satellite in the
%list PRNLIST.
%[satPositions, satClkCorr] = satpos(transmitTime, prnList, eph);
%
%   Inputs:
%       transmitTime  - transmission time: 1 by settings.numberOfChannels
%       prnList       - list of PRN-s to be processed
%       eph           - ephemeridies of satellites
%
%   Outputs:
%       satPositions  - positions of satellites (in ECEF system [X; Y; Z])
%       satClkCorr    - correction of satellites clocks in s

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

% GPS constatns

gpsPi          = 3.1415926535898;  % Pi used in the GPS coordinate 
                                   % system

%--- Constants for satellite position calculation -------------------------
Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM             = 3.986005e14;      % Earth's universal gravitational constant,
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
    dt = check_t(transmitTime(satNr) - eph(prn).t_oc);

    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - ...
                         eph(prn).T_GD;

    time = transmitTime(satNr) - satClkCorr(satNr);

%% Find satellite's position ----------------------------------------------

    % Restore semi-major axis
    a   = eph(prn).sqrtA * eph(prn).sqrtA;

    % Time correction
    tk  = check_t(time - eph(prn).t_oe);

    % Initial mean motion
    n0  = sqrt(GM / a^3);
    % Mean motion
    n   = n0 + eph(prn).deltan;

    % Mean anomaly
    M   = eph(prn).M_0 + n * tk;
    % Reduce mean anomaly to between 0 and 360 deg
    M   = rem(M + 2*gpsPi, 2*gpsPi);

    %Initial guess of eccentric anomaly
    E   = M;
    
    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end
    
    % Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*gpsPi, 2*gpsPi);

    % Relativistic correction
    dtr = F * eph(prn).e * eph(prn).sqrtA * sin(E);

    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);

    %Compute angle phi
    phi = nu + eph(prn).omega;
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*gpsPi);

    %Correct argument of latitude
    u = phi + ...
        eph(prn).C_uc * cos(2*phi) + ...
        eph(prn).C_us * sin(2*phi);		
    % Correct radius
    r = a * (1 - eph(prn).e*cos(E)) + ...
        eph(prn).C_rc * cos(2*phi) + ...
        eph(prn).C_rs * sin(2*phi);
    % Correct inclination
    i = eph(prn).i_0 + eph(prn).iDot * tk + ...
        eph(prn).C_ic * cos(2*phi) + ...
        eph(prn).C_is * sin(2*phi);
    
    % 2.9 SV position in orbital plane
    xk1 = cos(u)*r;
    yk1 = sin(u)*r;

    %Compute the angle between the ascending node and the Greenwich meridian
    Omega = eph(prn).omega_0 + (eph(prn).omegaDot - Omegae_dot)*tk - ...
            Omegae_dot * eph(prn).t_oe;
    %Reduce to between 0 and 360 deg
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

    %--- Compute satellite coordinates ------------------------------------
    xk = xk1 * cos(Omega) - yk1 * cos(i)*sin(Omega);
    yk = xk1 * sin(Omega) + yk1 * cos(i)*cos(Omega);
    zk = yk1 * sin(i);
    satPositions(1, satNr) = xk;
    satPositions(2, satNr) = yk;
    satPositions(3, satNr) = zk;

%% Include relativistic correction in clock correction --------------------
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + ...
                         eph(prn).a_f0 - ...
                         eph(prn).T_GD + dtr;       
                     
%% The following is to calculate sv velocity (currently not used in this version)
% %% Computation of SV velocity in ECEF -----------------------------------
%     % 2.12 Derivative of eccentric Anomaly ----------------------------
%     dE = n/(1-eph(prn).e *cos(E));
%     
%     % 2.13 Derivative of argument of Latitude --------------------------
%     dphi = sqrt(1 - eph(prn).e^2) * dE / (1-eph(prn).e *cos(E));
%     
%     % 2.14-2.15 Derivative of the following terms ------------------------   
%     % Derivative of argument of latitude
%     du = dphi + ...
%         2*dphi*(-eph(prn).C_uc * sin(2*phi) + ...
%         eph(prn).C_us * cos(2*phi));
%     
% 	% Derivative of radius
%     %Correct radius
%     dr = a * eph(prn).e * dE *sin(E) + ...
%         2*dphi*(-eph(prn).C_rc * sin(2*phi) + ...
%         eph(prn).C_rs * cos(2*phi));
%     
% 	% Derivative of inclination
%     %Correct inclination
%     di = eph(prn).iDot + ...
%         2*dphi*(-eph(prn).C_ic * sin(2*phi) + ...
%         eph(prn).C_is * cos(2*phi));
%     
%     % Derivative of Longitude of Ascending Node
%     dOmega = eph(prn).omegaDot - Omegae_dot;
%     
%     % 2.16 SV velocity in orbital plane ------------------------
%     dxk1 = dr*cos(u) - r*du*sin(u);
%     dyk1 = dr*sin(u) + r*du*cos(u);
%     
%     % 2.17 SV velocity in ECEF ------------------------
%     satVolocity(1, satNr) = -yk*dOmega - (dyk1*cos(i) - zk*di) * sin(Omega) + dxk1*cos(Omega);
%     satVolocity(2, satNr) = xk*dOmega  + (dyk1*cos(i) - zk*di) * cos(Omega) + dxk1*sin(Omega) ;
%     satVolocity(3, satNr) = dyk1*sin(i) + yk1*di*cos(i);
% 
% %% 4¡¢Include relativistic correction in clock rate correction  -------------
%     % Relativistic correction
%     dtrRat = F * eph(prn).e * eph(prn).sqrtA * cos(E) *dE;
%     
% 	% The clock drift is relative small, thus can be neglectde at most time.
%     satClkCorrRat(satNr) = 2* eph(prn).a_f2 * dt + eph(prn).a_f1 + dtrRat;                    
   
end % for satNr = 1 : numOfSatellites
