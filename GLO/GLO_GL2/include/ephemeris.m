function [eph, TOD] = ephemeris(navBiBinaryBitsSamples)
%Function decodes ephemerides and TOD from the given bit stream. The stream
%(array) in the parameter BITS must contain 1275 bits. The first element in
%the array must be the first bit of a string. The string ID of the
%first string in the array is not important.
%
%
%[eph, TOD] = ephemeris(bits)
%
%   Inputs:
%       bits        - bits of the navigation messages (15 strings).
%                   Type is character array and it must contain only
%                   characters '0' or '1'.
%   Outputs:
%       TOD         - Time Of Day (TOD) of the first string in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris
%
%--------------------------------------------------------------------------
%                           CU Multi-GNSS SDR  
%
% Copyright (C) Darius Plausinaitis and Kristin Larson
% Written by Darius Plausinaitis and Kristin Larson
%
% GLONASS modification by Jakob Almqvist
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
%$Id: ephemeris.m,v 1.1.2.7 2006/08/14 11:38:22 dpl Exp $

% This is in order to make sure variable 'eph' for each SV has a similar
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
eph = eph_structure_init();

%--- Group every 10 vales of bits into columns ------------------------
navBiBinaryBitsSamples = reshape(navBiBinaryBitsSamples, ...
    10, (size(navBiBinaryBitsSamples, 1) / 10));

%--- Sum all samples in the bits to get the best estimate -------------
navBiBinaryBits = sum(navBiBinaryBitsSamples);

%--- Now threshold and make 1 and 0 -----------------------------------
% The expression (navBits > 0) returns an array with elements set to 1
% if the condition is met and set to 0 if it is not met.
navBiBinaryBits = (navBiBinaryBits > 0);

%--- Convert from bi-binary to relative code -------------------------
relNavBits = (navBiBinaryBits(1:2:2549) - ...
    navBiBinaryBits(2:2:2550) +1 )./2;

%--- Convert from relative code to data sequence and checking bits ---
navBits(1) = 0;
navBits(2:1275) = xor(relNavBits(1:1274),relNavBits(2:1275));

%--- Convert from decimal to binary -----------------------------------
% The function ephemeris expects input in binary form. In Matlab it is
% a string array containing only "0" and "1" characters.
navBitsBin = dec2bin(navBits);

bits = navBitsBin(1:1275)';
%% Check if there is enough data ==========================================
if length(bits) < 1275
    error('The parameter BITS must contain 1275 bits!');
end

%% Check if the parameters are strings ====================================
if ~ischar(bits)
    error('The parameter BITS must be a character array!');
end

%% Decode all 15 strings ==================================================
for i = 1:15
    
    %--- "Cut" one string's bits ------------------------------------------
    string = bits(85*(i-1)+1 : 85*i);
    
    %--- Correct polarity of the data bits in the string ------------------
    string = checkPhase(string);
    
    %--- Decode the string id ---------------------------------------------
    % For more details on string contents please refer to GLONASS ICD.
    stringID = bin2dec(string(2:5));
    
    %--- Decode string based on the strings id ----------------------------
    % The task is to select the necessary bits and convert them to decimal
    % numbers. For more details on string contents please refer to GLONASS
    % ICD.
    switch stringID
        case 1  %--- It is string 1 ---------------------------------------
            % It contains TOD, flag P1 and x- coordinate, velocity and
            % acceleration.
            eph.P1     = bin2dec(string(8:9));
            if eph.P1 ~= 0
                eph.P1 = (eph.P1+1) * 15;                         %[flag]
            end
            eph.TOD    = bin2dec(string(10:14)) * 3600 + ...      %[sec]
                bin2dec(string(15:20)) * 60 + ...
                bin2dec(string(21)) * 30;
            eph.xDis   = (-1)^bin2dec(string(51)) * ...           %[km]
                bin2dec(string(52:77)) * 2^(-11);
            eph.xVel   = (-1)^bin2dec(string(22)) * ...           %[km/s]
                bin2dec(string(23:45)) * 2^(-20);
            eph.xAcc   = (-1)^bin2dec(string(46)) * ...           %[km/s^2]
                bin2dec(string(47:50)) * 2^(-30);
            
        case 2  %--- It is string 2 ---------------------------------------
            % It contains reference time for ephemeris, P2 flag, health
            % flag and y- coordinate, velocity and acceleration
            eph.P2     = bin2dec(string(9));                      %[flag]
            eph.B      = bin2dec(string(6));                      %[flag]
            eph.tb     = bin2dec(string(10:16)) * 15 * 60;        %[sec]
            eph.yDis   = (-1)^bin2dec(string(51)) * ...           %[km]
                bin2dec(string(52:77)) * 2^(-11);
            eph.yVel   = (-1)^bin2dec(string(22)) * ...           %[km/s]
                bin2dec(string(23:45)) * 2^(-20);
            eph.yAcc   = (-1)^bin2dec(string(46)) * ...           %[km/s^2]
                bin2dec(string(47:50)) * 2^(-30);
            
        case 3  %--- It is string 3 ---------------------------------------
            % It contains frequency offset, P flag, P3 flag, health flag
            % and z- coordinate, velocity and acceleration
            eph.P3     = bin2dec(string(6));                      %[flag]
            eph.gam  = (-1)^bin2dec(string(7)) * ...
                bin2dec(string(8:17)) * 2^(-40);
            eph.P      = bin2dec(string(19:20));                  %[flag]
            eph.health = bin2dec(string(21));
            eph.zDis   = (-1)^bin2dec(string(51)) * ...           %[km]
                bin2dec(string(52:77)) * 2^(-11);
            eph.zVel   = (-1)^bin2dec(string(22)) * ...           %[km/s]
                bin2dec(string(23:45)) * 2^(-20);
            eph.zAcc   = (-1)^bin2dec(string(46)) * ...           %[km/s^2]
                bin2dec(string(47:50)) * 2^(-30);
            
        case 4 %--- It is string 4 ----------------------------------------
            % It contains fourth part of ephemeris paramters
            eph.tau_n  = (-1)^bin2dec(string(6)) * ...            %[sec]
                bin2dec(string(7:27)) * 2^(-30);
            eph.dtau   = (-1)^bin2dec(string(28)) * ...           %[sec]
                bin2dec(string(29:32)) * 2^(-30);
            eph.E      = bin2dec(string(33:37));                  %[days]
            eph.P4     = bin2dec(string(52));                     %[flag]
            eph.FT     = bin2dec(string(53:56));
            eph.M      = bin2dec(string(76:77));
            eph.n      = bin2dec(string(71:75));
            eph.days   = bin2dec(string(60:70));              %[days]
            
        case 5 %--- It is string 5 ----------------------------------------
            % It constains year and GLONASS time scale correction to UTC.
            eph.N4     = bin2dec(string(50:54));
            eph.tau_c  = (-1)^bin2dec(string(17)) * ...
                bin2dec(string(18:48)) * 2^(-31);
            
        otherwise %--- All the other strings ------------------------------
            % Almanac
            % Not decoded at the moment.
            
    end % switch stringID ...
    
end % for all 15 srings ...

%% Compute the time of day (TOD) of the first string in the array ========
% Also correct the TOD. The transmitted TOD is the time reference for
% when the frame where the TOD was given (string 1) started. The TOD is
% re-calculated so it corresponds to the first string in the data block.

TOD = eph.TOD - (15-stringID) * 2;


