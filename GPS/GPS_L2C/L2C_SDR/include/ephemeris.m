function [eph, TOW] = ephemeris(navBitsBin,eph)
%Function decodes ephemerides and TOW from the given bit stream. The stream
%(array) in the parameter BITS must contain 1500 bits. The first element in
%the array must be the first bit of a subframe. The subframe ID of the
%first subframe in the array is not important.
%
%Function does not check parity!
%
%[eph, TOW] = ephemeris(bits, D30Star)
%
%   Inputs:
%       navBitsBin  - bits of the navigation messages.Type is character array 
%                   and it must contain only characters '0' or '1'.
%       eph         - The ephemeris for each PRN is decoded message by message.
%                   To prevent lost of previous decoded messages, the eph sturcture
%                   must be passed onto this function.
%   Outputs:
%       TOW         - Time Of Week (TOW) of the first sub-frame in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

%% Check if there is enough data ==========================================
if length(navBitsBin) < 300
    error('The parameter BITS must contain 1500 bits!');
end

%% Check if the parameters are strings ======================================
if ~ischar(navBitsBin)
    error('The parameter BITS must be a character array!');
end

% Pi used in the GPS coordinate system
gpsPi = 3.1415926535898; 

%% Decode messages needed to compute positionning ===========================

%--- Decode the message id ----------------------------------------------------------
% For more details on message contents please refer to IS-GPS200H.
messageID = bin2dec(navBitsBin(15:20));

%--- Decode messages based on the message id ----------------------------------------
% The task is to select the necessary bits and convert them to decimal
% numbers. For more details on message contents please refer to GPS
% ICD (IS-GPS200H.).
switch messageID
    % Message type 10 in conjunction with message type 11 provides users the
    % requisite data to calculate SV position.
    case 10  %--- It is Message Type 10 ---------------------------------------------
        % It contains first part of ephemeris parameters
        eph.idValid(1) = 10;
        % PRN
        eph.PRN  = bin2dec(navBitsBin(9:14));
        % Week No.
        eph.weekNumber  = bin2dec(navBitsBin(39:51));
        % Top Time of Ephemeris prediction
        eph.T_op       = bin2dec(navBitsBin(55:65))*300;
        % L2 health
        eph.health      = bin2dec(navBitsBin(53));
        % ED Accuracy Index
        eph.URA_ED      = twosComp2dec(navBitsBin(66:70));
        % Ephemeris data reference time of week
        eph.t_oe        = bin2dec(navBitsBin(71:81)) * 300;  
        % Semi-major axis difference at reference time
        eph.deltaA      = twosComp2dec(navBitsBin(82:107)) * 2^(-9) ;
        % Change rate in semi-major axis
        eph.ADot        = twosComp2dec(navBitsBin(108:132)) * 2^(-21);
        % Mean Motion difference from computed value at reference time
        eph.delta_n_0   = twosComp2dec(navBitsBin(133:149)) * 2^(-44)* gpsPi;
        % Rate of mean motion difference from computed value
        eph.delta_n_0Dot= twosComp2dec(navBitsBin(150:172)) * 2^(-57)* gpsPi;
        % Mean anomaly at reference time
        eph.M_0         = twosComp2dec(navBitsBin(173:205)) * 2^(-32) * gpsPi;
        % Eccentricity
        eph.e           = bin2dec(navBitsBin(206:238))* 2^(-34);
        % Argument of perigee
        eph.omega       = twosComp2dec(navBitsBin(239:271))* 2^(-32) * gpsPi;
    
    case 11  %--- It is Message Type 11 ---------------------------------------------
        % It contains second part of ephemeris parameters
        eph.idValid(2)  = 11;
        % PRN
        eph.PRN  = bin2dec(navBitsBin(9:14));
        % Ephemeris data reference time of week
        eph.t_oe        = bin2dec(navBitsBin(39:49)) * 300;
        % Longitude of Ascending Node of Orbit Plane at Weekly Epoch
        eph.omega_0     = twosComp2dec(navBitsBin(50:82))* 2^(-32) * gpsPi;
        % Inclination angle at reference time
        eph.i_0         = twosComp2dec(navBitsBin(83:115))* 2^(-32) * gpsPi;
        % Rate of right ascension difference
        eph.delta_omegaDot  = twosComp2dec(navBitsBin(116:132)) * 2^(-44) * gpsPi;
        % Rate of inclination angle
        eph.i_0Dot      = twosComp2dec(navBitsBin(133:147)) * 2^(-44) * gpsPi;       
        % Amplitude of the sine harmonic correction term to the angle of inclination
        eph.C_is        = twosComp2dec(navBitsBin(148:163)) * 2^(-30);
        % Amplitude of the cosine harmonic correction term to the angle of inclination
        eph.C_ic        = twosComp2dec(navBitsBin(164:179)) * 2^(-30);
        % Amplitude of the sine correction term to the orbit radius
        eph.C_rs        = twosComp2dec(navBitsBin(180: 203)) * 2^(-8);
        % Amplitude of the cosine correction term to the orbit radius
        eph.C_rc        = twosComp2dec(navBitsBin(204:227)) * 2^(-8);
        % Amplitude of the sine harmonic correction term to the argument of latitude
        eph.C_us        = twosComp2dec(navBitsBin(228:248)) * 2^(-30);
        % Amplitude of the cosine harmonic correction term to the argument of latitude
        eph.C_uc        = twosComp2dec(navBitsBin(249:269)) * 2^(-30);
                                
    case 30 %--- It is Message Type 30 ----------------------------------------------
        % It contains Clock, IONO & Group Delay
        eph.idValid(3) = 30;
        % PRN
        eph.PRN  = bin2dec(navBitsBin(9:14));
        % Clock Data Reference Time of Week
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        % SV Clock Bias Correction Coefficient
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        % SV Clock Drift Correction Coefficient
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        % SV Clock Drift Rate Correction Coefficient
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % The group delay differential correction terms
        eph.T_GD        = twosComp2dec(navBitsBin(128:140)) * 2^(-35);
        eph.ISC_L2C     = twosComp2dec(navBitsBin(154:166)) * 2^(-35);
        % The ionospheric parameters
        eph.alpha0      = twosComp2dec(navBitsBin(193:200)) * 2^(-30);
        eph.alpha1      = twosComp2dec(navBitsBin(201:208)) * 2^(-27);
        eph.alpha2      = twosComp2dec(navBitsBin(209:216)) * 2^(-24);
        eph.alpha3      = twosComp2dec(navBitsBin(217:224)) * 2^(-24);
        eph.beta0       = twosComp2dec(navBitsBin(225:232)) * 2^(11);
        eph.beta1       = twosComp2dec(navBitsBin(233:240)) * 2^(14);
        eph.beta2       = twosComp2dec(navBitsBin(241:248)) * 2^(16);
        eph.beta3       = twosComp2dec(navBitsBin(249:256)) * 2^(16);
        
    case 31 %--- It is Message Type 31 ----------------------------------------------
        % It contains Clock & Reduced Almanac
        eph.idValid(4) = 31;
        eph.PRN  = bin2dec(navBitsBin(9:14));
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % Other terms not decoded at the moment...
        
    case 32 %--- It is Message Type 32 ----------------------------------------------
        % It contains Clock & EOP
        eph.idValid(5) = 32;
        eph.PRN  = bin2dec(navBitsBin(9:14));
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % Other terms not decoded at the moment...
        
    case 33 %--- It is Message Type 33 ----------------------------------------------
        % It contains Clock & UTC
        eph.idValid(6) = 33;
        eph.PRN  = bin2dec(navBitsBin(9:14));
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % Other terms not decoded at the moment...
        
    case 34 %--- It is Message Type 34 ----------------------------------------------
        % It contains Clock & Differential Correction
        eph.idValid(7)  = 34;
        eph.PRN  = bin2dec(navBitsBin(9:14));
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % Other terms not decoded at the moment...
        
    case 35 %--- It is Message Type 35 ----------------------------------------------
        % It contains Clock & GGTO
        eph.idValid(8)  = 35;
        eph.PRN  = bin2dec(navBitsBin(9:14));
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % Other terms not decoded at the moment...
        
    case 36 %--- It is Message Type 36 ----------------------------------------------
        % It contains Clock & Text
        eph.idValid(9) = 36;
        eph.PRN  = bin2dec(navBitsBin(9:14));
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % Other terms not decoded at the moment...
        
    case 37 %--- It is Message Type 37 ----------------------------------------------
        % It contains Clock & Midi Almanac
        eph.idValid(10) = 37;
        eph.PRN  = bin2dec(navBitsBin(9:14));
        eph.t_oc        = bin2dec(navBitsBin(61:71)) * 300;
        eph.a_f0        = twosComp2dec(navBitsBin(72:97)) * 2^(-35);
        eph.a_f1        = twosComp2dec(navBitsBin(98:117)) * 2^(-48);
        eph.a_f2        = twosComp2dec(navBitsBin(118:127)) * 2^(-60);
        % Other terms not decoded at the moment...
        
    otherwise % Other message types include: ----------------------------------------      
        % Mainly Reduced & Midi Almanac,UTC parameters and so on
        % Not decoded at the moment.
        eph.idValid(11) = messageID;
        % PRN
        eph.PRN  = bin2dec(navBitsBin(9:14));
       
end % switch subframeID ...

%% Compute the time of week (TOW) of the first message ======================
% Also correct the TOW. The transmitted TOW is actual TOW of the next
% message and we need the TOW of the start of current message in this data block.
% It is the same as HOW in NAV.
TOW = bin2dec(navBitsBin(21:37)) * 6 - 12;