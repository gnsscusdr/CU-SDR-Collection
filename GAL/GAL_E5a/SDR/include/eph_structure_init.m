function eph = eph_structure_init ()
% Initialize the Ephemeris data before decoding.===========================

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------

%--- SVID , Clock correction, SISA, Ionospheric correction,
% BGD, GST, Signal health and Data validity status ------------------------
eph.SVID      = [];       % 6 bits
eph.IODnav1   = [];       % 10 bits
eph.t_oc      = [];       % 14 bits in s
eph.a_f0      = [];       % 31 bit in s
eph.a_f1      = [];       % 21 bit in s/s
eph.a_f2      = [];       % 6 bit in s/(s^2)
eph.a_i0      = [];       % 11 bit
eph.a_i1      = [];       % 11 bit
eph.a_i2      = [];       % 14 bit
eph.iono_SF1  = [];       % 1 bit
eph.iono_SF2  = [];       % 1 bit
eph.iono_SF3  = [];       % 1 bit
eph.iono_SF4  = [];       % 1 bit
eph.iono_SF5  = [];       % 1 bit
eph.BGD_E1E5a = [];       % 10 bit
eph.E5a_HS    = [];       % 2 bit  E5a signal Health Status (0-OK)
eph.WN        = [];       % 12 bit
eph.TOW       = [];       % 20 bit
eph.E5a_DVS   = [];       % 1 bit E5a Data validity status (0-valid)

%--- Ephemeris (1/3) and GST -------------------------------------
eph.IODnav2   = [];       % 10 bit
eph.M_0       = [];       % 32 bit in rad
eph.OmegaDot  = [];       % 24 bit in rad
eph.e         = [];       % 32 bit
eph.sqrtA     = [];       % 32 bit in m^0.2
eph.Omega_0   = [];       % 32 bit in rad
eph.iDot      = [];       % 14 bit in rad

%--- Ephemeris (2/3) and GST -------------------------------
eph.IODnav3   = [];       % 10 bit
eph.i_0       = [];       % 32 bit in rad
eph.omega     = [];       % 32 bit in rad
eph.deltan    = [];       % 16 bit in rad
eph.CUC       = [];       % 16 bit in rad
eph.CUS       = [];       % 16 bit in rad
eph.CRC       = [];       % 16 bit in m
eph.CRS       = [];       % 16 bit in m
eph.t_oe      = [];       % 14 bit in s

%--- Ephemeris (3/3), GST-UTC conversion, GST-GPS conversion and TO W
eph.IODnav4   = [];       % 10 bit
eph.CIC       = [];       % 16 bit in rad
eph.CIS       = [];       % 16 bit in rad
eph.A0        = [];       % 32 bit
eph.A1        = [];       % 24 bit
eph.delt_LS   = [];       % 8 bit
eph.t_ot      = [];       % 8 bit
eph.WN_ot     = [];       % 8 bit
eph.WN_LSF    = [];       % 8 bit
eph.DN        = [];       % 3 bit
eph.delt_LSF  = [];       % 8 bit
eph.t_og      = [];       % 8 bit
eph.A0_G      = [];       % 16 bit
eph.A1_G      = [];       % 12 bit
eph.WN_og     = [];       % 6 bit

% ephemeris decoded flag --------------------------------------------------
eph.flag         = [];
