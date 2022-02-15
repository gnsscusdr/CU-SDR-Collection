function [eph] = ephemeris(navBitsBin,eph)
%Function decodes ephemerides and TOW from the given bit stream. The stream
%(array) in the parameter BITS must contain 288 bits. The first element in
%the array must be the first bit of a subframe.
%
%
%[eph] = ephemeris(navBitsBin,eph)
%
%   Inputs:
%       navBitsBin  - bits of the navigation messages.Type is character array
%                   and it must contain only characters '0' or '1'.
%       eph         - The ephemeris for each PRN is decoded message by message.
%                   To prevent lost of previous decoded messages, the eph sturcture
%                   must be passed onto this function.
%   Outputs:
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

%% Preparation for date message decoding ============================
if length(navBitsBin) < 288
    error('The parameter BITS must contain 1500 bits!');
end

% Check if the parameters are strings
if ~ischar(navBitsBin)
    error('The parameter BITS must be a character array!');
end

% 'bits' should be row vector for 'bin2dec' function.
[a, b] = size(navBitsBin);
if a > b
    navBitsBin = navBitsBin';
end

% Pi used in the GPS coordinate system
gpsPi = 3.1415926535898;

% Decode PRN
PRN = bin2dec(navBitsBin(1:6));
if (PRN < 1) || (PRN >63)
    return
end

% Decode the message id
MesType = bin2dec(navBitsBin(7:12));

%%  Decode messages based on the message id =========================
% The task is to select the necessary bits and convert them to decimal
% numbers. For more details on message contents please refer to BDS-3
% ICD (BDS-SIS-ICD-B2a-1.0).
switch MesType
    % Message type 10 in conjunction with message type 11 provides users
    % the requisite data to calculate SV position.
    case 10  %--- It is Message Type 10 -----------------------------------
        % It contains first part of ephemeris parameters
        eph.idValid(1) = 10;
        % PRN
        eph.PRN  = PRN;
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        % Week No.
        eph.WN  = bin2dec(navBitsBin(31:43));
        % DIF
        eph.DIF  = bin2dec(navBitsBin(44));
        % SIF
        eph.SIF  = bin2dec(navBitsBin(45));
        % AIF
        eph.AIF  = bin2dec(navBitsBin(46));
        % Ephemeris data reference time of week
        eph.t_oe        = bin2dec(navBitsBin(62:72)) * 300;
        % Satellite type
        SatType     = bin2dec(navBitsBin(73:74));
        if (SatType == 1)
            eph.SatType = 'GEO';
        elseif (SatType == 2)
            eph.SatType = 'IGSO';
        elseif (SatType == 3)
            eph.SatType = 'MEO';
        end
        % Semi-major axis difference at reference time
        eph.deltaA      = twosComp2dec(navBitsBin(75:100)) * 2^(-9) ;
        % Change rate in semi-major axis
        eph.ADot        = twosComp2dec(navBitsBin(101:125)) * 2^(-21);
        % Mean Motion difference from computed value at reference time
        eph.delta_n_0   = twosComp2dec(navBitsBin(126:142)) * 2^(-44)* gpsPi;
        % IRate of mean motion difference from computed value
        eph.delta_n_0Dot= twosComp2dec(navBitsBin(143:165)) * 2^(-57)* gpsPi;
        % Mean anomaly at reference time
        eph.M_0         = twosComp2dec(navBitsBin(166:198)) * 2^(-32) * gpsPi;
        % Eccentricity
        eph.e           = bin2dec(navBitsBin(199:231))* 2^(-34);
        % Argument of perigee
        eph.omega       = twosComp2dec(navBitsBin(232:264))* 2^(-32) * gpsPi;
        
    case 11  %--- It is Message Type 11 -----------------------------------
        % It contains second part of ephemeris parameters
        eph.idValid(2)  = 11;
        % PRN
        eph.PRN  = PRN;
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        % HS
        eph.HS  = bin2dec(navBitsBin(31:32));
        % DIF
        eph.DIF  = bin2dec(navBitsBin(33));
        % SIF
        eph.SIF  = bin2dec(navBitsBin(34));
        % AIF
        eph.AIF  = bin2dec(navBitsBin(36));
        % Longitude of Ascending Node of Orbit Plane at Weekly Epoch
        eph.omega_0     = twosComp2dec(navBitsBin(43:75))* 2^(-32) * gpsPi;
        % Inclination angle at reference time
        eph.i_0         = twosComp2dec(navBitsBin(76:108))* 2^(-32) * gpsPi;
        % Rate of right ascension difference
        eph.omegaDot  = twosComp2dec(navBitsBin(109:127)) * 2^(-44) * gpsPi;
        % Rate of inclination angle
        eph.i_0Dot      = twosComp2dec(navBitsBin(128:142)) * 2^(-44) * gpsPi;
        % Amplitude of the sine harmonic correction term to the angle of inclination
        eph.C_is        = twosComp2dec(navBitsBin(143:158)) * 2^(-30);
        % Amplitude of the cosine harmonic correction term to the angle of inclination
        eph.C_ic        = twosComp2dec(navBitsBin(159:174)) * 2^(-30);
        % Amplitude of the sine correction term to the orbit radius
        eph.C_rs        = twosComp2dec(navBitsBin(175: 198)) * 2^(-8);
        % Amplitude of the cosine correction term to the orbit radius
        eph.C_rc        = twosComp2dec(navBitsBin(199:222)) * 2^(-8);
        % Amplitude of the sine harmonic correction term to the argument of latitude
        eph.C_us        = twosComp2dec(navBitsBin(223:243)) * 2^(-30);
        % Amplitude of the cosine harmonic correction term to the argument of latitude
        eph.C_uc        = twosComp2dec(navBitsBin(244:264)) * 2^(-30);
        
    case 30 %--- It is Message Type 30 ------------------------------------
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        % It contains Clock, IONO & Group Delay
        eph.idValid(3) = 30;
        % PRN
        eph.PRN  = PRN;
        % Clock Data Reference Time of Week
        eph.t_oc        = bin2dec(navBitsBin(43:53)) * 300;
        % SV Clock Bias Correction Coefficient
        eph.a_0        = twosComp2dec(navBitsBin(54:78)) * 2^(-34);
        % SV Clock Drift Correction Coefficient
        eph.a_1        = twosComp2dec(navBitsBin(79:100)) * 2^(-50);
        % SV Clock Drift Rate Correction Coefficient
        eph.a_2        = twosComp2dec(navBitsBin(101:111)) * 2^(-66);
        % MSB 2 bits
        eph.IODC_MSB2   = bin2dec(navBitsBin(112:113));
        % LSB 8 bits
        eph.IODC_LSB8   = bin2dec(navBitsBin(114:121));
        
        % Group delay differential of the B2a pilot componet
        eph.T_GDB2ap        = twosComp2dec(navBitsBin(122:133)) * 2^(-34);
        % Group delay differential between the B1C data and pilot components
        eph.ISC_B2ad        = twosComp2dec(navBitsBin(134:145)) * 2^(-34);

        % The ionospheric parameters
        eph.alpha1      = bin2dec(navBitsBin(146:155)) * 2^(-3);
        eph.alpha2      = twosComp2dec(navBitsBin(156:163)) * 2^(-3);
        eph.alpha3      = bin2dec(navBitsBin(164:171)) * 2^(-3);
        eph.alpha4      = bin2dec(navBitsBin(172:179)) * 2^(-3);
        eph.alpha5       = bin2dec(navBitsBin(180:187)) * 2^(-3);
        eph.alpha6       = twosComp2dec(navBitsBin(188:195)) * 2^(-3);
        eph.alpha7       = twosComp2dec(navBitsBin(196:203)) * 2^(-3);
        eph.alpha8       = twosComp2dec(navBitsBin(204:211)) * 2^(-3);
        eph.alpha9       = twosComp2dec(navBitsBin(212:219)) * 2^(-3);
        
        % Group delay differential of the B1C pilot component
        eph.T_GDB1Cp        = twosComp2dec(navBitsBin(220:231)) * 2^(-34);
        
    case 31 %--- It is Message Type 31 ------------------------------------
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        % It contains Clock & Reduced Almanac
        eph.idValid(4) = 31;
        % ORN NO.
        eph.PRN  = PRN;
        % Clock Data Reference Time of Week
        eph.t_oc        = bin2dec(navBitsBin(43:53)) * 300;
        % SV Clock Bias Correction Coefficient
        eph.a_0        = twosComp2dec(navBitsBin(54:78)) * 2^(-34);
        % SV Clock Drift Correction Coefficient
        eph.a_1        = twosComp2dec(navBitsBin(79:100)) * 2^(-50);
        % SV Clock Drift Rate Correction Coefficient
        eph.a_2        = twosComp2dec(navBitsBin(101:111)) * 2^(-66);
        % MSB 2 bits
        eph.IODC_MSB2   = bin2dec(navBitsBin(112:113));
        % LSB 8 bits
        eph.IODC_LSB8   = bin2dec(navBitsBin(114:121));
        % Other terms not decoded at the moment...
        
    case 32 %--- It is Message Type 32 ------------------------------------
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        % It contains Clock & EOP
        eph.idValid(5) = 32;
        % PRN NO.
        eph.PRN  = PRN;
        % Clock Data Reference Time of Week
        eph.t_oc        = bin2dec(navBitsBin(43:53)) * 300;
        % SV Clock Bias Correction Coefficient
        eph.a_0        = twosComp2dec(navBitsBin(54:78)) * 2^(-34);
        % SV Clock Drift Correction Coefficient
        eph.a_1        = twosComp2dec(navBitsBin(79:100)) * 2^(-50);
        % SV Clock Drift Rate Correction Coefficient
        eph.a_2        = twosComp2dec(navBitsBin(101:111)) * 2^(-66);
        % MSB 2 bits
        eph.IODC_MSB2   = bin2dec(navBitsBin(112:113));
        % LSB 8 bits
        eph.IODC_LSB8   = bin2dec(navBitsBin(114:121));
        % Other terms not decoded at the moment...
        
    case 33 %--- It is Message Type 33 ------------------------------------
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        % It contains Clock & UTC
        eph.idValid(6) = 33;
        % PRN No.
        eph.PRN  = PRN;
        % Clock Data Reference Time of Week
        eph.t_oc        = bin2dec(navBitsBin(43:53)) * 300;
        % SV Clock Bias Correction Coefficient
        eph.a_0        = twosComp2dec(navBitsBin(54:78)) * 2^(-34);
        % SV Clock Drift Correction Coefficient
        eph.a_1        = twosComp2dec(navBitsBin(79:100)) * 2^(-50);
        % SV Clock Drift Rate Correction Coefficient
        eph.a_2        = twosComp2dec(navBitsBin(101:111)) * 2^(-66);
        % MSB 2 bits
        eph.IODC_MSB2   = bin2dec(navBitsBin(112:113));
        % LSB 8 bits
        eph.IODC_LSB8   = bin2dec(navBitsBin(114:121));
        % BGTO --------------------------------------
        % GNSS ID
        eph.GNSS_ID   = bin2dec(navBitsBin(112:114));
        eph.WN_0BGTO   = bin2dec(navBitsBin(115:127));
        eph.t_0BGTO   = bin2dec(navBitsBin(128:143))* 2^(4);
        eph.A_0BGTO   = twosComp2dec(navBitsBin(144:159))* 2^(-35);
        eph.A_1BGTO   = twosComp2dec(navBitsBin(160:172))* 2^(-51);
        eph.A_2BGTO   = twosComp2dec(navBitsBin(173:179))* 2^(-68);
        % Other terms not decoded at the moment...
        
    case 34 %--- It is Message Type 34 ------------------------------------
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        % It contains Clock & Differential Correction
        eph.idValid(7)  = 34;
        % PRN No.
        eph.PRN  = PRN;
        % Clock Data Reference Time of Week
        eph.t_oc        = bin2dec(navBitsBin(65:75)) * 300;
        % SV Clock Bias Correction Coefficient
        eph.a_0        = twosComp2dec(navBitsBin(76:100)) * 2^(-34);
        % SV Clock Drift Correction Coefficient
        eph.a_1        = twosComp2dec(navBitsBin(101:122)) * 2^(-50);
        % SV Clock Drift Rate Correction Coefficient
        eph.a_2        = twosComp2dec(navBitsBin(123:133)) * 2^(-66);
        % MSB 2 bits
        eph.IODC_MSB2   = bin2dec(navBitsBin(134:135));
        % LSB 8 bits
        eph.IODC_LSB8   = bin2dec(navBitsBin(136:143));
        % BDT-UTC ------------------------------------------
        eph.A_0UTC        = twosComp2dec(navBitsBin(123:133)) * 2^(-35);
        eph.A_1UTC        = twosComp2dec(navBitsBin(123:133)) * 2^(-51);
        eph.A_2UTC        = twosComp2dec(navBitsBin(123:133)) * 2^(-68);
        eph.delta_t_LS    = twosComp2dec(navBitsBin(123:133));
        eph.t_ot          = bin2dec(navBitsBin(123:133)) * 2^(4);
        eph.WN_ot         = bin2dec(navBitsBin(123:133));
        eph.WN_LSF        = bin2dec(navBitsBin(123:133));
        eph.DN            = bin2dec(navBitsBin(123:133));
        eph.delta_t_LSF   = twosComp2dec(navBitsBin(123:133));
        % Other terms not decoded at the moment...
        
    otherwise % Other message types include: ------------------------------
        % Mainly Reduced & Midi Almanac,UTC parameters and so on
        % Not decoded at the moment.
        eph.idValid(8) = MesType;
        % PRN
        eph.PRN  = PRN;
        % SOW
        if isempty(eph.SOW)
            eph.SOW  = bin2dec(navBitsBin(13:30)) * 3;
        end
        
end % switch MesType ...
