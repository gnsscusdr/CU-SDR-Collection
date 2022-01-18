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
eph.idValid(1:5) = zeros(1,5);
% PRN
eph.PRN  = [];

%--- It is subframe 1 -------------------------------------
% It contains WN, SV clock corrections, health and accuracy
eph.weekNumber  = [];
eph.accuracy    = [];
eph.health      = [];
eph.T_GD        = [];
eph.IODC        = [];
eph.t_oc        = [];
eph.a_f2        = [];
eph.a_f1        = [];
eph.a_f0        = [];

%--- It is subframe 2 -------------------------------------
% It contains first part of ephemeris parameters
eph.IODE_sf2    = [];
eph.C_rs        = [];
eph.deltan      = [];
eph.M_0         = [];
eph.C_uc        = [];
eph.e           = [];
eph.C_us        = [];
eph.sqrtA       = [];
eph.t_oe        = [];

%--- It is subframe 3 -------------------------------------
% It contains second part of ephemeris parameters
eph.C_ic        = [];
eph.omega_0     = [];
eph.C_is        = [];
eph.i_0         = [];
eph.C_rc        = [];
eph.omega       = [];
eph.omegaDot    = [];
eph.IODE_sf3    = [];
eph.iDot        = [];

%--- It is subframe 4 -------------------------------------
% Almanac, ionospheric model, UTC parameters.
% SV health (PRN: 25-32).
% Not decoded at the moment.

%--- It is subframe 5 -------------------------------------
% SV almanac and health (PRN: 1-24).
% Almanac reference week number and time.
% Not decoded at the moment.

% Tow of first decoded subframe
eph.TOW         = [];