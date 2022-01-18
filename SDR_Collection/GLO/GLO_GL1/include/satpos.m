function [satPositions, satClkCorr] = satpos(transmitTime, satList, ...
    eph, tau_c) 
%SATPOS Calculation of X,Y,Z satellites coordinates at TRANSMITTIME for
%given ephemeris EPH. Coordinates are calculated for each satellite in the
%list SATLIST.
%[satPositions, satClkCorr] = satpos(transmitTime, satList, eph, settings);
%
%   Inputs:
%       transmitTime  - transmission time
%       satList       - list of satellites to be processed
%       eph           - ephemeridies of satellites
%       tau_c         - difference between GLONASS time and UTC(SU) time
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
%GLONASS modification by Jakob Almqvist
%
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $

%% Initialize constants ===================================================
numOfSatellites = size(satList, 2);

%--- Constants for satellite position calculation -------------------------
omega          = 7.292115e-5;      % Earth rotation rate, [rad/s]
my             = 3.986004418e14;   % Earth's universal
                                   % gravitational parameter,
                                   % [m^3/s^2]
a              = 6.378136e6;       % Semi-major axis of Earth, [m]
J02            = 1.082657e-3;      % Second zonal harmonic of the
                                   % geopotential
%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for satNr = 1 : numOfSatellites
    
    chn = satList(satNr);
    
%% Find initial satellite clock correction --------------------------------

    %--- Find time difference ---------------------------------------------
    dt = check_t(transmitTime(chn) - eph(chn).tb);
    
    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = -(eph(chn).tau_n + tau_c - eph(chn).gam * dt);
    
    %--- Find integration time --------------------------------------------
    time    = dt - satClkCorr(satNr);   
    
    %--- Check if to integrate forward or backward ------------------------
    if time < 0
        tau      = -60;
    else
        tau      = 60;
    end
    
    
%% Find satellite's position ----------------------------------------------

    % x,y,z- coordinates
    x0  = eph(chn).xDis * 1e3;
    y0  = eph(chn).yDis * 1e3;
    z0  = eph(chn).zDis * 1e3;

    % x,y,z- velocities
    Vx0 = eph(chn).xVel * 1e3;
    Vy0 = eph(chn).yVel * 1e3;
    Vz0 = eph(chn).zVel * 1e3;
    
    % x,y,z- accelerations
    Ax0 = eph(chn).xAcc * 1e3;
    Ay0 = eph(chn).yAcc * 1e3;
    Az0 = eph(chn).zAcc * 1e3;
    
    %--- Start integrating using Runge-Kutta -----------------------------
    for numberOfIntegrations = tau:tau:time+tau
    
        % Check if last integration step. If last integration step, make 
        % one more step that has the remaining time length...
        if abs(numberOfIntegrations) > abs(time)
            
            % ...if there is more time left to integrate...
            if mod(time,tau) ~= 0
                tau = mod(time,tau);
                
            % ...otherwise make a zero step.
            else
                tau = 0;
            end
        end
        
       
        r0  = sqrt(x0^2+y0^2+z0^2);
        
        %--- D1 -----------------------------------------------------------
        dVx_1  = - my / r0^3 * x0 - 3/2 * J02 * my * a^2 / r0^5 * x0 * ...
              (1 - 5 * z0^2 / r0^2) + omega^2 * x0 + 2 * omega * Vy0 + Ax0;
        dVy_1  = - my / r0^3 * y0 - 3/2 * J02 * my * a^2 / r0^5 * y0 * ...
              (1 - 5 * z0^2 / r0^2) + omega^2 * y0 - 2 * omega * Vx0 + Ay0;
        dVz_1  = - my / r0^3 * z0 - 3/2 * J02 * my * a^2 / r0^5 * z0 * ...
              (3 - 5 * z0^2 / r0^2) + Az0;
        
        dx_1  = Vx0;
        dy_1  = Vy0;
        dz_1  = Vz0;  
          
        %------ Xn  + 1/2 * tau * D1 --------------------------------------  
        Vx_1  = Vx0 + 1/2 * tau * dVx_1;
        Vy_1  = Vy0 + 1/2 * tau * dVy_1;
        Vz_1  = Vz0 + 1/2 * tau * dVz_1;
               
        x_1   = x0  + 1/2 * tau * dx_1;
        y_1   = y0  + 1/2 * tau * dy_1;
        z_1   = z0  + 1/2 * tau * dz_1;
        
        r_1  = sqrt( x_1^2 + y_1^2 + z_1^2 );
        
        %--- D2 -----------------------------------------------------------
        dVx_2  = - my / r_1^3 * x_1 - 3/2 * J02 * my * a^2 / r_1^5 * x_1 * ...
              (1 - 5 * z_1^2 / r_1^2) + omega^2 * x_1 + 2 * omega * Vy_1 + Ax0;
        dVy_2  = - my / r_1^3 * y_1 - 3/2 * J02 * my * a^2 / r_1^5 * y_1 * ...
              (1 - 5 * z_1^2 / r_1^2) + omega^2 * y_1 - 2 * omega * Vx_1 + Ay0;
        dVz_2  = - my / r_1^3 * z_1 - 3/2 * J02 * my * a^2 / r_1^5 * z_1 * ...
              (3 - 5 * z_1^2 / r_1^2) + Az0;
        
        dx_2  = Vx_1;
        dy_2  = Vy_1;
        dz_2  = Vz_1;  
        
        %------ Xn  + 1/2 * tau * D2 --------------------------------------  
        Vx_2  = Vx0 + 1/2 * tau * dVx_2;
        Vy_2  = Vy0 + 1/2 * tau * dVy_2;
        Vz_2  = Vz0 + 1/2 * tau * dVz_2;
               
        x_2   = x0  + 1/2 * tau * dx_2;
        y_2   = y0  + 1/2 * tau * dy_2;
        z_2   = z0  + 1/2 * tau * dz_2;
        
        r_2  = sqrt( x_2^2 + y_2^2 + z_2^2 );
        
        %--- D3 -----------------------------------------------------------
        dVx_3  = - my / r_2^3 * x_2 - 3/2 * J02 * my * a^2 / r_2^5 * x_2 * ...
              (1 - 5 * z_2^2 / r_2^2) + omega^2 * x_2 + 2 * omega * Vy_2 + Ax0;
        dVy_3  = - my / r_2^3 * y_2 - 3/2 * J02 * my * a^2 / r_2^5 * y_2 * ...
              (1 - 5 * z_2^2 / r_2^2) + omega^2 * y_2 - 2 * omega * Vx_2 + Ay0;
        dVz_3  = - my / r_2^3 * z_2 - 3/2 * J02 * my * a^2 / r_2^5 * z_2 * ...
              (3 - 5 * z_2^2 / r_2^2) + Az0;
        
        dx_3  = Vx_2;
        dy_3  = Vy_2;
        dz_3  = Vz_2;  
        
        %------ Xn  + tau * D3 -------------------------------------------  
        Vx_3  = Vx0 + tau * dVx_3;
        Vy_3  = Vy0 + tau * dVy_3;
        Vz_3  = Vz0 + tau * dVz_3;
               
        x_3   = x0  + tau * dx_3;
        y_3   = y0  + tau * dy_3;
        z_3   = z0  + tau * dz_3;
        
        r_3  = sqrt( x_3^2 + y_3^2 + z_3^2 );
        
        %--- D4 -----------------------------------------------------------
        dVx_4  = - my / r_3^3 * x_3 - 3/2 * J02 * my * a^2 / r_3^5 * x_3 * ...
              (1 - 5 * z_3^2 / r_3^2) + omega^2 * x_3 + 2 * omega * Vy_3 + Ax0;
        dVy_4  = - my / r_3^3 * y_3 - 3/2 * J02 * my * a^2 / r_3^5 * y_3 * ...
              (1 - 5 * z_3^2 / r_3^2) + omega^2 * y_3 - 2 * omega * Vx_3 + Ay0;
        dVz_4  = - my / r_3^3 * z_3 - 3/2 * J02 * my * a^2 / r_3^5 * z_3 * ...
              (3 - 5 * z_3^2 / r_3^2) + Az0;
        
        dx_4  = Vx_3;
        dy_4  = Vy_3;
        dz_4  = Vz_3;  
        
        %------------------------------------------------------------------
        
        Vx0  = Vx0 + 1/6 * tau * ( dVx_1 + 2 * dVx_2 + ...
            2 * dVx_3 + dVx_4 );
        Vy0  = Vy0 + 1/6 * tau * ( dVy_1 + 2 * dVy_2 + ...
            2 * dVy_3 + dVy_4 );
        Vz0  = Vz0 + 1/6 * tau * ( dVz_1 + 2 * dVz_2 + ...
            2 * dVz_3 + dVz_4 );
             
        x0   = x0 + 1/6 * tau * ( dx_1 + 2 * dx_2 + ...
            2 * dx_3 + dx_4 );
        y0   = y0 + 1/6 * tau * ( dy_1 + 2 * dy_2 + ...
            2 * dy_3 + dy_4 );
        z0   = z0 + 1/6 * tau * ( dz_1 + 2 * dz_2 + ...
            2 * dz_3 + dz_4 );
        
    end
    
    %--- Compute satellite coordinates -------------------------------------
    satPositions(1, satNr) = x0;
    satPositions(2, satNr) = y0;
    satPositions(3, satNr) = z0;

    %--- Transform from PZ-90 to WGS84 -------------------------------------  
     satPositions(:,satNr) = [1       -1.6e-6  0;
                              1.6e-6   1       0;
                              0        0       1] * satPositions(:,satNr);

end % for satNr = 1 : numOfSatellites