function eph= eph_structure_init()
% This is in order to make sure variable 'eph' for each SV has a similar
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (c) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

%--- It is string 1 ---------------------------------------
% It contains TOD, flag P1 and x- coordinate, velocity and
% acceleration.
eph.P1     = [];      %[flag]
eph.TOD    = [];      %[sec]
eph.xDis   = [];      %[km]
eph.xVel   = [];      %[km/s]
eph.xAcc   = [];      %[km/s^2]

%--- It is string 2 ---------------------------------------
% It contains reference time for ephemeris, P2 flag, health
% flag and y- coordinate, velocity and acceleration
eph.P2     = [];      %[flag]
eph.B      = [];      %[flag]
eph.tb     = [];      %[sec]
eph.yDis   = [];      %[km]
eph.yVel   = [];      %[km/s]
eph.yAcc   = [];      %[km/s^2]

%--- It is string 3 ---------------------------------------
% It contains frequency offset, P flag, P3 flag, health flag
% and z- coordinate, velocity and acceleration
eph.P3     = [];      %[flag]
eph.gam  = [];
eph.P      = [];      %[flag]
eph.health = [];
eph.zDis   = [];      %[km]
eph.zVel   = [];      %[km/s]
eph.zAcc   = [];      %[km/s^2]

%--- It is string 4 ----------------------------------------
% It contains fourth part of ephemeris paramters
eph.tau_n  = [];      %[sec]
eph.dtau   = [];      %[sec]
eph.E      = [];      %[days]
eph.P4     = [];      %[flag]
eph.FT     = [];
eph.M      = [];
eph.n      = [];
eph.days   = [];      %[days]

%--- It is string 5 ----------------------------------------
% It constains year and GLONASS time scale correction to UTC.
eph.N4     = [];
eph.tau_c  = [];
eph.validFlag = [];