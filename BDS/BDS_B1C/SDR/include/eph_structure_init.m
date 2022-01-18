function eph= eph_structure_init()
% This is in order to make sure variable 'eph' for each SV has a similar
% structure when only one or even none of the three requisite sub-frames
% is decoded for a given PRN.
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
%$Id: eph_structure_init.m,v 1.1.2.7 2018/03/06 11:38:22 dpl Exp $

%% ===== Decode the first subframe ==================================
eph.PRN  = [];
% SOH
eph.SOH  = [];

%% ===== Decode the second subframe =================================
% Week No.
eph.WN  = [];
% HOW
eph.HOW  = [];
% IODC
eph.IODC = [];
% IODE
eph.IODE = [];

% --- Ephemeris I ---------------------------------------------------------
% Ephemeris data reference time of week
eph.t_oe        = [];
% Satellite Type
eph.SatType     = [];
% Semi-major axis difference at reference time
eph.deltaA      = [];
% Change rate in semi-major axis
eph.ADot        = [];
% Mean Motion difference from computed value at reference time
eph.delta_n_0   = [];
% IRate of mean motion difference from computed value
eph.delta_n_0Dot= [];
% Mean anomaly at reference time
eph.M_0         = [];
% Eccentricity
eph.e           = [];
% Argument of perigee
eph.omega       = [];

% --- Ephemeris II --------------------------------------------------------
% Longitude of Ascending Node of Orbit Plane at Weekly Epoch
eph.omega_0     = [];
% Inclination angle at reference time
eph.i_0         = [];
% Rate of right ascension difference
eph.omegaDot  = [];
% Rate of inclination angle
eph.i_0Dot      = [];
% Amplitude of the sine harmonic correction term to the angle of inclination
eph.C_is        = [];
% Amplitude of the cosine harmonic correction term to the angle of inclination
eph.C_ic        = [];
% Amplitude of the sine correction term to the orbit radius
eph.C_rs        = [];
% Amplitude of the cosine correction term to the orbit radius
eph.C_rc        = [];
% Amplitude of the sine harmonic correction term to the argument of latitude
eph.C_us        = [];
% Amplitude of the cosine harmonic correction term to the argument of latitude
eph.C_uc        = [];

% --- SV clock error parameters -------------------------------------------
% Clock Data Reference Time of Week
eph.t_oc        = [];
% SV Clock Bias Correction Coefficient
eph.a_0        = [];
% SV Clock Drift Correction Coefficient
eph.a_1        = [];
% SV Clock Drift Rate Correction Coefficient
eph.a_2        = [];

% --- Remaining parts of the second subframe-------------------------------
eph.T_GDB2ap        = [];
eph.ISC_B1Cd        = [];
eph.T_GDB1Cp        = [];


%% ===== Decode the third subframe ==================================
% PageID
eph.PageID1   = [];
eph.PageID3   = [];
% Heath state
eph.HS       = [];
% DIF
eph.DIF      = [];
% SIF
eph.SIF      = [];
% AIF
eph.AIF      = [];
% SISMAI
eph.SISMAI   = [];

% other parts of subframe 3

% The ionospheric parameters ----------------------------------------------
eph.alpha1      = [];
eph.alpha2      = [];
eph.alpha3      = [];
eph.alpha4      = [];
eph.alpha5       = [];
eph.alpha6       = [];
eph.alpha7       = [];
eph.alpha8       = [];
eph.alpha9       = [];

% BDT-UTC -----------------------------------------------------------------
eph.A_0UTC        = [];
eph.A_1UTC        = [];
eph.A_2UTC        = [];
eph.delta_t_LS    = [];
eph.t_ot          = [];
eph.WN_ot         = [];
eph.WN_LSF        = [];
eph.DN            = [];
eph.delta_t_LSF   = [];

% BGTO --------------------------------------------------------------------
% GNSS ID
eph.GNSS_ID   = [];
eph.WN_0BGTO   = [];
eph.t_0BGTO   = [];
eph.A_0BGTO   = [];
eph.A_1BGTO   = [];
eph.A_2BGTO   = [];

% All required B-CNAV1 sub-frame message has been decoded
eph.flag = [];
% Time of Week
eph.TOW = [];

