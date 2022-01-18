function eph= eph_structure_init()
% This is in order to make sure variable 'eph' for each SV has a similar 
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
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
%$Id: ephemeris.m,v 1.1.2.7 2017/03/06 11:38:22 dpl Exp $

% Flags for message data decoding. 0 indicates decoding fail, 1 is successful
% decoding. idValid(1:2) for message type 10 and 11, idValid(3:10) for message
% type 30:37, idValid(11) for others.
eph.idValid(1:11) = zeros(1,11);
% PRN
eph.PRN  = [];
%% Message type 10 ==========================================================
% Week No.
eph.weekNumber  = [];
% Top Time of Ephemeris prediction
eph.T_op       = [];
% L2 health
eph.health      = [];
% ED Accuracy Index
eph.URA_ED      = [];
% Ephemeris data reference time of week
eph.t_oe        = [];
% Semi-major axis difference at reference time
eph.deltaA      = [];
% Change rate in semi-major axis
eph.ADot        = [];
% Mean Motion difference from computed value at reference time
eph.delta_n_0   = [];
% Rate of mean motion difference from computed value
eph.delta_n_0Dot= [];
% Mean anomaly at reference time
eph.M_0         = [];
% Eccentricity
eph.e           = [];
% Argument of perigee
eph.omega       = [];
%% Message type 11 ==========================================================
% Longitude of Ascending Node of Orbit Plane at Weekly Epoch
eph.omega_0     = [];
% Inclination angle at reference time
eph.i_0         = [];
% Rate of right ascension difference
eph.delta_omegaDot  = [];
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
%% Message type 30 ==========================================================
% Clock Data Reference Time of Week
eph.t_oc        = [];
% SV Clock Bias Correction Coefficient
eph.a_f0        = [];
% SV Clock Drift Correction Coefficient
eph.a_f1        = [];
% SV Clock Drift Rate Correction Coefficient
eph.a_f2        = [];
% The group delay differential correction terms
eph.T_GD        = [];
eph.ISC_L5I     = [];
% The ionospheric parameters
eph.alpha0      = [];
eph.alpha1      = [];
eph.alpha2      = [];
eph.alpha3      = [];
eph.beta0       = [];
eph.beta1       = [];
eph.beta2       = [];
eph.beta3       = [];
% Tow of first decoded subframe
eph.TOW         = [];