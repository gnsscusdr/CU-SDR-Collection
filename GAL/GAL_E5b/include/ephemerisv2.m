function [eph, validWord] = ephemerisv2(navWord, wordType, eph)
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------

%Indicator the requisite messages are all decoded

%% Define required constants ==============================================
% Pi used in the Galileo coordinate system (same as for GPS)
galPi = 3.1415926535898;
validWord = 0;
%% Decode Messages ========================================================

%--- Decode sub-frame based on the message type -----------------------
switch wordType
    case 1   %--- Ephemeris (1/4) -------------------------------------
        eph.IODnav   = bin2dec(navWord(7:16));
        eph.t_oe      = bin2dec(navWord(17:30))        * 60;    % in s
        eph.M_0       = twosComp2dec(navWord(31:62))   * 2^(-31) * galPi; % in rad
        eph.e         = bin2dec(navWord(63:94))        * 2^(-33);
        eph.sqrtA     = bin2dec(navWord(95:126))       * 2^(-19);  % in m^0.2
        validWord = 1;
        
    case 2   %--- Ephemeris (2/4) -------------------------------------
        eph.IODnav   = bin2dec(navWord(7:16));
        eph.Omega_0   = twosComp2dec(navWord(17:48))   * 2^(-31) * galPi; % in rad
        eph.i_0       = twosComp2dec(navWord(49:80))   * 2^(-31) * galPi; % in rad
        eph.omega     = twosComp2dec(navWord(81:112))  * 2^(-31) * galPi; % in rad
        eph.iDot      = twosComp2dec(navWord(113:126)) * 2^(-43) * galPi; % in rad
        validWord = 1;
        
    case 3   %--- Ephemeris (3/4), SISA -------------------------------
        eph.IODnav   = bin2dec(navWord(7:16));
        eph.OmegaDot  = twosComp2dec(navWord(17:40))   * 2^(-43) * galPi; % in rad
        eph.deltan    = twosComp2dec(navWord(41:56))   * 2^(-43) * galPi; % in rad
        eph.CUC       = twosComp2dec(navWord(57:72))   * 2^(-29);         % in rad
        eph.CUS       = twosComp2dec(navWord(73:88))   * 2^(-29);         % in rad
        eph.CRC       = twosComp2dec(navWord(89:104))  * 2^(-5);          % in m
        eph.CRS       = twosComp2dec(navWord(105:120)) * 2^(-5);          % in m
        %eph.SISA      = navWord(121:128) % content not currently defined
        validWord = 1;
        
    case 4   %--- SVID, Ephemeris (4/4), Clock Correction -------------
        eph.IODnav   = bin2dec(navWord(7:16));
        eph.SVID      = bin2dec(navWord(17:22));
        eph.CIC       = twosComp2dec(navWord(23:38))   * 2^(-29); % in rad
        eph.CIS       = twosComp2dec(navWord(39:54))   * 2^(-29); % in rad
        eph.t_oc      = bin2dec(navWord(55:68))        * 60;      % in s
        eph.a_f0      = twosComp2dec(navWord(69:99))   * 2^(-34); % in s
        eph.a_f1      = twosComp2dec(navWord(100:120)) * 2^(-46); % in s/s
        eph.a_f2      = twosComp2dec(navWord(121:126)) * 2^(-59); % in s/(s^2)
        %             eph.a_f2      = 0; % --------------------------------------------------------------------------
        validWord = 1;
        
    case 5   %--- Iono, BGD, Signal Health, Data Validity, GST --------
        eph.a_i0      = bin2dec(navWord(7:17))         * 2^(-2);
        eph.a_i1      = twosComp2dec(navWord(18:28))   * 2^(-8);
        eph.a_i2      = twosComp2dec(navWord(29:42))   * 2^(-15);
        eph.iono_SF1  = bin2dec(navWord(43));
        eph.iono_SF2  = bin2dec(navWord(44));
        eph.iono_SF3  = bin2dec(navWord(45));
        eph.iono_SF4  = bin2dec(navWord(46));
        eph.iono_SF5  = bin2dec(navWord(47));
        eph.BGD_E1E5a = twosComp2dec(navWord(48:57))   * 2^(-32);
        eph.BGD_E1E5b = twosComp2dec(navWord(58:67))   * 2^(-32);
        eph.E5b_HS    = bin2dec(navWord(68:69));  % E5b signal Health Status (0-OK)
        eph.E1b_HS    = bin2dec(navWord(70:71));  % E1-B/C signal Health Status (0-OK)
        eph.E5b_DVS   = bin2dec(navWord(72));     % E5b Data validity status (0-valid)
        eph.E1B_DVS   = bin2dec(navWord(73));     % E1-B Data validity status (0-valid)
        eph.WN        = bin2dec(navWord(74:85));
        eph.TOW       = bin2dec(navWord(86:105));
        % Correct TOW to time for first page part
        validWord = 1;
        
    case 6   %--- GST-UTC Conversion ----------------------------------
        eph.A0        = twosComp2dec(navWord(7:38))    * 2^(-30);
        eph.A1        = twosComp2dec(navWord(39:62))   * 2^(-50);
        eph.delt_LS   = twosComp2dec(navWord(63:70));
        eph.t_ot      = bin2dec(navWord(71:78))        * 3600;
        eph.WN_ot     = bin2dec(navWord(79:86));
        eph.WN_LSF    = bin2dec(navWord(87:94));
        eph.DN        = bin2dec(navWord(95:97));
        eph.delt_LSF  = twosComp2dec(navWord(98:105));
        validWord = 1;
        
        %case 7   %--- Almanac -------------------------------------
        
        %case 8   %--- Almanac -------------------------------------
        
        %case 9   %--- Almanac -------------------------------------
        
    case 10  %--- GPS-GST Conversion ---------------------------------
        eph.A0_G      = twosComp2dec(navWord(87:102))  * 2^(-35);
        eph.A1_G      = twosComp2dec(navWord(103:114)) * 2^(-51);
        eph.t_og      = bin2dec(navWord(115:122))      * 3600;
        eph.WN_og     = bin2dec(navWord(123:128));
        validWord = 1;
end % switch word type

end % all pages
