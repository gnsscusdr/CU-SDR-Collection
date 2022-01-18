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
A_REF          = 26559710;         % Reference value for semi-major axis, [meter]
omegaDot_REF   = -2.6e-9 * gpsPi;  % Reference value for rate of right ascension,
                                   % [rad/s]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);
satVolocity  = zeros(3, numOfSatellites);
satClkCorrRat = zeros(1, numOfSatellites);

%% Process each satellite =================================================
for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    
    %% 1 Transmitting time correction --------------------------------------
    
    % ---- Find time difference -----------------------------------------
    % Accounting for beginning or end of week crossover
    dt = check_t(transmitTime(satNr) - eph(prn).t_oc);
    
    % ---- Calculate clock correction ----------------------------------
    % Relativistic correction can be only computed thereafter with E and
    % dE, so it is omitted and not included here.
    satClkCorr(satNr) = (eph(prn).a_f2 * dt + eph(prn).a_f1) * dt + eph(prn).a_f0;
    
    % ---- Include group delay for clock correction --------------------
    % If message type 30 including 'T_GD' and 'ISC_L2C' is decoded, then
    % make group delay correction. If the IF signal is not long enough
    % for decoding message type 30, then this correction will be omitted
    % for PVT calculaion. For more details on message contents please refer
    % to GPS ICD (IS-GPS200H) p171.
    if(eph(prn).idValid(3) == 30)
        satClkCorr(satNr) = satClkCorr(satNr) - ...
            eph(prn).T_GD + eph(prn).ISC_L5I;
    end
     
    % ---- Time from clock reference time ------------------------------
    time = transmitTime(satNr) - satClkCorr(satNr);
    
    %% 2 Compute satellite's position --------------------------------------
    % This algorithm can be found in "Xie Gang, Principles of GNSS: GPS, 
    % Glonass and Galileo" P248
    
    % ---- 2.1 Time from ephemeris reference time ----------------------
    % Accounting for beginning or end of week crossover
    tk  = check_t(time - eph(prn).t_oe);
    
    % ---- 2.2 Compute semi-major axis ---------------------------------
    % Semi-Major Axis at reference time
    % ** This is different with from NAV in L1 C/A **
    A_0   = A_REF + eph(prn).deltaA;
    % Semi-Major Axis
    A     = A_0 + eph(prn).ADot * tk;
    
    % ---- 2.3 Compute mean motion ------------------------------------
    %Initial mean motion
    n0  = sqrt(GM / A_0^3);
    % Mean motion difference from computed value
    % ** This is different with from NAV in L1 C/A **
    delta_n = eph(prn).delta_n_0 + 0.5 * eph(prn).delta_n_0Dot *tk;
    % Mean motion
    n   = n0 + delta_n;
    
    % ---- 2.4 Mean Anomaly computation --------------------------------
    %Mean anomaly
    M   = eph(prn).M_0 + n * tk;
    %Reduce mean anomaly to between 0 and 360 deg
    M   = rem(M + 2*gpsPi, 2*gpsPi);
    
    % ---- 2.5 Eccentric Anomaly --------------------------------------
    %Initial guess of eccentric anomaly
    E   = M;
    %Iteratively compute eccentric anomaly
    for ii = 1:10
        E_old   = E;
        E       = M + eph(prn).e * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);
        
        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end
    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*gpsPi, 2*gpsPi);
    
    % ---- 2.6 True Anomaly  ------------------------------------------
    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(prn).e^2) * sin(E), cos(E)-eph(prn).e);
    
    % ---- 2.7 Argument of Latitude -----------------------------------
    %Compute angle phi
    phi = nu + eph(prn).omega;
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*gpsPi);
    
    % ---- 2.8-2.9 Argument of latitude, radius and inclination -------
    %Correct argument of latitude
    u = phi + eph(prn).C_uc * cos(2*phi) + eph(prn).C_us * sin(2*phi);
    %Correct radius
    r = A * (1 - eph(prn).e*cos(E)) + eph(prn).C_rc * cos(2*phi) + ...
        eph(prn).C_rs * sin(2*phi);
    %Correct inclination
    i = eph(prn).i_0 + eph(prn).i_0Dot * tk + ...
        eph(prn).C_ic * cos(2*phi) + eph(prn).C_is * sin(2*phi);
    
    % ---- 2.10 SV position in orbital plane --------------------------
    xk1 = cos(u)*r;
    yk1 = sin(u)*r;
    
    % ---- 2.11 Longitude of Ascending Node ----------------------------
    % Rate of right ascension
    % ** This is different with from NAV in L1 C/A **
    omegaDot = omegaDot_REF + eph(prn).delta_omegaDot;
    % Corrected Longitude of Ascending Node
    Omega = eph(prn).omega_0 + (omegaDot - Omegae_dot)*tk - ...
        Omegae_dot * eph(prn).t_oe;
    % Reduce to between 0 and 360 deg
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);
    
    % ---- 2.12 Compute satellite coordinates in ECEF -----------------
    xk = xk1 * cos(Omega) - yk1 * cos(i)*sin(Omega);
    yk = xk1 * sin(Omega) + yk1 * cos(i)*cos(Omega);
    zk = yk1 * sin(i);
    
    satPositions(1, satNr) = xk;
    satPositions(2, satNr) = yk;
    satPositions(3, satNr) = zk;
    
    %% 3 Include relativistic correction in clock correction ----
    % Relativistic correction term
    dtr = F * eph(prn).e * sqrt(A) * sin(E);
    satClkCorr(satNr) = satClkCorr(satNr) + dtr;
    
%% The following is to calculate sv velocity (currently not used in this version)
%     %% 4 Computation of SV velocity in ECEF ---------------------
%     % This algorithm can be found in "Xie Gang, Principles of GNSS: GPS, 
%     % Glonass and Galileo" P249
%     % ---- 4.1-4.2 Derivative of eccentric Anomaly ---------------------
%     dE = n/(1-eph(prn).e *cos(E));
%     
%     % ---- 4.3-4.4 Derivative of argument of Latitude ------------------
%     dphi = sqrt(1 - eph(prn).e^2) * dE / (1-eph(prn).e *cos(E));
%     
%     % ---- 4.5-4.6 Derivative of the following terms  ------------------
%     % Derivative of argument of latitude
%     du = dphi + 2*dphi*(-eph(prn).C_uc * sin(2*phi) + ...
%                                       eph(prn).C_us * cos(2*phi));
%     
%     % Derivative of radius
%     dr = A * eph(prn).e * dE *sin(E) + 2*dphi*(-eph(prn).C_rc * sin(2*phi) + ...
%                                       eph(prn).C_rs * cos(2*phi));
%     
%     % Derivative of inclination
%     di = eph(prn).i_0Dot + 2*dphi*(-eph(prn).C_ic * sin(2*phi) + ...
%                                        eph(prn).C_is * cos(2*phi));
%     
%     % ---- 4.7 SV velocity in orbital plane ---------------------------
%     dxk1 = dr*cos(u) - r*du*sin(u);
%     dyk1 = dr*sin(u) + r*du*cos(u);
%     
%     % ---- 4.8 Derivative of Longitude of Ascending Node
%     dOmega = omegaDot - Omegae_dot;
% 
%     % ---- 4.9 SV velocity in ECEF ------------------------------------
%     satVolocity(1, satNr) = -yk*dOmega - (dyk1*cos(i) - ...
%                               zk*di) * sin(Omega) + dxk1*cos(Omega);
%     satVolocity(2, satNr) = xk*dOmega  + (dyk1*cos(i) - ...
%                               zk*di) * cos(Omega) + dxk1*sin(Omega) ;
%     satVolocity(3, satNr) = dyk1*sin(i) + yk1*di*cos(i);
%     
%     %% 5 Include relativistic correction in clock rate correction  ----------
%     % Relativistic correction
%     dtrRat = F * eph(prn).e * sqrt(A) * cos(E) *dE;
%     
%     % The clock drift is relative small, thus can be neglectde at most time.
%     satClkCorrRat(satNr) = 2* eph(prn).a_f2 * dt + eph(prn).a_f1 + dtrRat;
    
end % for satNr = 1 : numOfSatellites
