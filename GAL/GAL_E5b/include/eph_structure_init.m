function eph = eph_structure_init ()
% Initialize the Ephemeris data before decoding.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------
eph = struct();
eph.TOW = [];
%--- Ephemeris (1/4) ------------------------------------------------------
eph.IODnav   = [];
eph.t_oe      = [];
eph.M_0       = [];
eph.e         = [];
eph.sqrtA     = [];

%--- Ephemeris (2/4) ------------------------------------------------------
eph.Omega_0   = [];
eph.i_0       = [];
eph.omega     = [];
eph.iDot      = [];

%--- Ephemeris (3/4), SISA ------------------------------------------------
eph.OmegaDot  = [];
eph.deltan    = [];
eph.CUC       = [];
eph.CUS       = [];
eph.CRC       = [];
eph.CRS       = [];
%eph.SISA      = []; % content not currently defined

%--- SVID, Ephemeris (4/4), Clock Correction ------------------------------
eph.SVID      = [];
eph.CIC       = [];
eph.CIS       = [];
eph.t_oc      = [];
eph.a_f0      = [];
eph.a_f1      = [];
eph.a_f2      = [];

%--- Iono, BGD, Signal Health, Data Validity, GST -------------------------
eph.a_i0      = [];
eph.a_i1      = [];
eph.a_i2      = [];
eph.iono_SF1  = [];
eph.iono_SF2  = [];
eph.iono_SF3  = [];
eph.iono_SF4  = [];
eph.iono_SF5  = [];
eph.BGD_E1E5a = [];
eph.BGD_E1E5b = [];
eph.E5b_HS    = [];
eph.E1b_HS    = [];
eph.E5b_DVS   = [];
eph.E1B_DVS   = [];
eph.WN        = [];
% Correct TOW to time for first page part
eph.TOW = [];

%--- GST-UTC Conversion ---------------------------------------------------
eph.A0        = [];
eph.A1        = [];
eph.delt_LS   = [];
eph.t_ot      = [];
eph.WN_ot     = [];
eph.WN_LSF    = [];
eph.DN        = [];
eph.delt_LSF  = [];

%--- GPS-GST Conversion ---------------------------------------------------
eph.A0_G      = [];
eph.A1_G      = [];
eph.t_og      = [];
eph.WN_og     = [];

% ephemeris decoded flag --------------------------------------------------
eph.flag         = 0;
