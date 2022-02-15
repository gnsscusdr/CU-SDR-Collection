function [eph] = ephemeris(navBitsBin,eph)
%Function decodes ephemerides and TOW from the given bit stream. The stream
%(array) in the parameter BITS must contain 1500 bits. The first element in
%the array must be the first bit of a subframe. The subframe ID of the
%first subframe in the array is not important.
%
%Function does not check parity!
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
%       TOW         - Time Of Week (TOW) of the first sub-frame in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos

% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and 
% implementation of an open-source BDS-3 B1C/B2a SDR receiver. 
% GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
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

% For more details on message contents please refer to BDS-SIS-ICD-B1C-1.0.


%% Check if the parameters are strings ==============================
if ~ischar(navBitsBin)
    error('The parameter BITS must be a character array!');
end

% 'bits' should be row vector for 'bin2dec' function.
[a, b] = size(navBitsBin);
if a > b
    navBitsBin = navBitsBin';
end

% Pi used in the GPS coordinate system
bdsPi = 3.1415926535898;

%% ===== Decode the first subframe ==================================
eph.PRN = bin2dec(navBitsBin(1:6));
if (eph.PRN < 1) || (eph.PRN > 63)
    return
end

if isempty(eph.flag)
    % Seconds Of Hour (SOH)
    eph.SOH  = bin2dec(navBitsBin(7:14)) * 18;
    
    %% ===== Decode the second subframe =================================
    subFra2Bit = navBitsBin(15:53);
    % Week Number (WN)
    eph.WN  = bin2dec(subFra2Bit(1:13));
    % Hours Of Week
    eph.HOW  = bin2dec(subFra2Bit(14:21));
    % IODC
    eph.IODC = bin2dec(subFra2Bit(22:32));
    % IODE
    eph.IODE = bin2dec(subFra2Bit(32:39));
    
    % --- Ephemeris I ---------------------------------------------------------
    subFra2Bit = navBitsBin(54:256);
    % Ephemeris data reference time of week
    eph.t_oe        = bin2dec(subFra2Bit(1:11)) * 300;
    % Satellite Type
    SatType     = bin2dec(subFra2Bit(12:13));
    if (SatType == 1)
        eph.SatType = 'GEO';
    elseif (SatType == 2)
        eph.SatType = 'IGSO';
    elseif (SatType == 3)
        eph.SatType = 'MEO';
    end
    % Semi-major axis difference at reference time
    eph.deltaA      = twosComp2dec(subFra2Bit(14:39)) * 2^(-9) ;
    % Change rate in semi-major axis
    eph.ADot        = twosComp2dec(subFra2Bit(40:64)) * 2^(-21);
    % Mean Motion difference from computed value at reference time
    eph.delta_n_0   = twosComp2dec(subFra2Bit(65:81)) * 2^(-44)* bdsPi;
    % IRate of mean motion difference from computed value
    eph.delta_n_0Dot= twosComp2dec(subFra2Bit(82:104)) * 2^(-57)* bdsPi;
    % Mean anomaly at reference time
    eph.M_0         = twosComp2dec(subFra2Bit(105:137)) * 2^(-32) * bdsPi;
    % Eccentricity
    eph.e           = bin2dec(subFra2Bit(138:170))* 2^(-34);
    % Argument of perigee
    eph.omega       = twosComp2dec(subFra2Bit(171:203))* 2^(-32) * bdsPi;
    
    % --- Ephemeris II --------------------------------------------------------
    subFra2Bit = navBitsBin(257:478);
    % Longitude of Ascending Node of Orbit Plane at Weekly Epoch
    eph.omega_0     = twosComp2dec(subFra2Bit(1:33))* 2^(-32) * bdsPi;
    % Inclination angle at reference time
    eph.i_0         = twosComp2dec(subFra2Bit(34:66))* 2^(-32) * bdsPi;
    % Rate of right ascension difference
    eph.omegaDot  = twosComp2dec(subFra2Bit(67:85)) * 2^(-44) * bdsPi;
    % Rate of inclination angle
    eph.i_0Dot      = twosComp2dec(subFra2Bit(86:100)) * 2^(-44) * bdsPi;
    % Amplitude of the sine harmonic correction term to the angle
    % of inclination
    eph.C_is        = twosComp2dec(subFra2Bit(101:116)) * 2^(-30);
    % Amplitude of the cosine harmonic correction term to the angle
    % of inclination
    eph.C_ic        = twosComp2dec(subFra2Bit(117:132)) * 2^(-30);
    % Amplitude of the sine correction term to the orbit radius
    eph.C_rs        = twosComp2dec(subFra2Bit(133: 156)) * 2^(-8);
    % Amplitude of the cosine correction term to the orbit radius
    eph.C_rc        = twosComp2dec(subFra2Bit(157:180)) * 2^(-8);
    % Amplitude of the sine harmonic correction term to the argument
    % of latitude
    eph.C_us        = twosComp2dec(subFra2Bit(181:201)) * 2^(-30);
    % Amplitude of the cosine harmonic correction term to the argument
    % of latitude
    eph.C_uc        = twosComp2dec(subFra2Bit(202:222)) * 2^(-30);
    
    % --- SV clock error parameters -------------------------------------------
    subFra2Bit = navBitsBin(479:547);
    % Clock Data Reference Time of Week
    eph.t_oc        = bin2dec(subFra2Bit(1:11)) * 300;
    % SV Clock Bias Correction Coefficient
    eph.a_0        = twosComp2dec(subFra2Bit(12:36)) * 2^(-34);
    % SV Clock Drift Correction Coefficient
    eph.a_1        = twosComp2dec(subFra2Bit(37:58)) * 2^(-50);
    % SV Clock Drift Rate Correction Coefficient
    eph.a_2        = twosComp2dec(subFra2Bit(59:69)) * 2^(-66);
    
    % --- Remaining parts of the second subframe-------------------------------
    subFra2Bit = navBitsBin(548:583);
    % Group delay differential of the B2a pilot compone
    eph.T_GDB2ap        = twosComp2dec(subFra2Bit(1:12)) * 2^(-34);
    % Group delay differential between the B1C data and pilot components
    eph.ISC_B1Cd        = twosComp2dec(subFra2Bit(13:24)) * 2^(-34);
    % Group delay differential of the B1C pilot component
    eph.T_GDB1Cp        = twosComp2dec(subFra2Bit(25:36)) * 2^(-34);
end

%% ===== Decode the third subframe ==================================
subFra2Bit = navBitsBin(615:878);
% PageID
PageID   = bin2dec(subFra2Bit(1:6));
% Heath state
eph.HS       = bin2dec(subFra2Bit(7:8));
% data integrity flag (DIF)
eph.DIF      = bin2dec(subFra2Bit(9));
% signal integrity flag (SIF)
eph.SIF      = bin2dec(subFra2Bit(10));
% accuracy integrity flag (AIF)
eph.AIF      = bin2dec(subFra2Bit(11));
% space monitoring accuracy index (SISMAI)
eph.SISMAI   = bin2dec(subFra2Bit(12:15));

% other parts of subframe 3
if PageID == 1
    eph.PageID1 = 1;
    % The ionospheric parameters ------------------------------------------
    tempBit = subFra2Bit(43:116);
    eph.alpha1      = bin2dec(tempBit(1:10)) * 2^(-3);
    eph.alpha2      = twosComp2dec(tempBit(11:18)) * 2^(-3);
    eph.alpha3      = bin2dec(tempBit(19:26)) * 2^(-3);
    eph.alpha4      = bin2dec(tempBit(27:34)) * 2^(-3);
    eph.alpha5       = bin2dec(tempBit(35:42)) * 2^(-3);
    eph.alpha6       = twosComp2dec(tempBit(43:50)) * 2^(-3);
    eph.alpha7       = twosComp2dec(tempBit(51:58)) * 2^(-3);
    eph.alpha8       = twosComp2dec(tempBit(59:66)) * 2^(-3);
    eph.alpha9       = twosComp2dec(tempBit(67:74)) * 2^(-3);
    
    % BDT-UTC -------------------------------------------------------------
    tempBit = subFra2Bit(117:213);
    % Bias coefficient of BDT time scale relative to UTC time scale
    eph.A_0UTC        = twosComp2dec(tempBit(1:16)) * 2^(-35);
    % Drift coefficient of BDT time scale relative to UTC time scale
    eph.A_1UTC        = twosComp2dec(tempBit(17:29)) * 2^(-51);
    % Drift rate coefficient of BDT time scale relative to UTC time scale
    eph.A_2UTC        = twosComp2dec(tempBit(30:36)) * 2^(-68);
    % Current or past leap second count
    eph.delta_t_LS    = twosComp2dec(tempBit(37:44));
    % Reference time of week
    eph.t_ot          = bin2dec(tempBit(45:60)) * 2^(4);
    % Reference week number
    eph.WN_ot         = bin2dec(tempBit(61:73));
    % Leap second reference week number
    eph.WN_LSF        = bin2dec(tempBit(74:86));
    % Leap second reference day number
    eph.DN            = bin2dec(tempBit(87:89));
    % Current or future leap second count
    eph.delta_t_LSF   = twosComp2dec(tempBit(90:97));
    
elseif PageID == 3
    eph.PageID3 = 3;
    % BDT-GNSS Time Offset (BGTO) -----------------------------------------
    tempBit = subFra2Bit(159:226);
    % GNSS type identification
    eph.GNSS_ID   = bin2dec(tempBit(1:3));
    % Reference week number
    eph.WN_0BGTO   = bin2dec(tempBit(4:16));
    % Reference time of week
    eph.t_0BGTO   = bin2dec(tempBit(17:32))* 2^(4);
    % Bias coefficient of BDT time scale relative to GNSS time scale
    eph.A_0BGTO   = twosComp2dec(tempBit(33:48))* 2^(-35);
    % Drift coefficient of BDT time scale relative to GNSS time scale
    eph.A_1BGTO   = twosComp2dec(tempBit(49:61))* 2^(-51);
    % Drift rate coefficient of BDT time scale relative to GNSS time scale
    eph.A_2BGTO   = twosComp2dec(tempBit(62:68))* 2^(-68);
end

if isempty(eph.flag)
    eph.TOW = eph.HOW * 3600 + eph.SOH;
end

% All required NAV data has been decoded
eph.flag = 1;

