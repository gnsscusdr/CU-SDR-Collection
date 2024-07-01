function [eph, validWord, wordType] = ephemeris(navWord, eph)

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------

%--- Convert navigation data word to binary ---------------------------

%--- Decode the message type ------------------------------------------
wordType = bin2dec(navWord(1:6));

% Galileo Pi definition
galPi = 3.1415926535898;
%
validWord = 0;
%--- Decode sub-frame based on the message type -----------------------
switch wordType
    case 1   %--- SVID , Clock correction, SISA, Ionospheric correction,
        % BGD, GST, Signal health and Data validity status
        eph.SVID      = bin2dec(navWord(7:12));                     % 6 bits
        eph.IODnav1   = bin2dec(navWord(13:22));                    % 10 bits
        eph.t_oc      = bin2dec(navWord(23:36))        * 60;        % 14 bits in s
        eph.a_f0      = twosComp2dec(navWord(37:67))   * 2^(-34);   % 31 bit in s
        eph.a_f1      = twosComp2dec(navWord(68:88)) * 2^(-46);     % 21 bit in s/s
        eph.a_f2      = twosComp2dec(navWord(89:94)) * 2^(-59);     % 6 bit in s/(s^2)
        %             eph.SISA_E1_E5a = bin2dec(navWord(95:102));                 % 8 bit           eph.SISA      = navWord(121:128) % content not currently defined
        eph.a_i0      = bin2dec(navWord(103:113))         * 2^(-2);    % 11 bit
        eph.a_i1      = twosComp2dec(navWord(114:124))   * 2^(-8);    % 11 bit
        eph.a_i2      = twosComp2dec(navWord(125:138))   * 2^(-15);   % 14 bit
        eph.iono_SF1  = bin2dec(navWord(139));      % 1 bit
        eph.iono_SF2  = bin2dec(navWord(140));      % 1 bit
        eph.iono_SF3  = bin2dec(navWord(141));      % 1 bit
        eph.iono_SF4  = bin2dec(navWord(142));      % 1 bit
        eph.iono_SF5  = bin2dec(navWord(143));      % 1 bit
        eph.BGD_E1E5a = twosComp2dec(navWord(144:153))   * 2^(-32);  % 10 bit
        eph.E5a_HS    = bin2dec(navWord(154:155));   % 2 bit  E5a signal Health Status (0-OK)
        eph.WN        = bin2dec(navWord(156:167));   % 12 bit
        eph.TOW       = bin2dec(navWord(168:187));  % 20 bit
        eph.E5a_DVS   = bin2dec(navWord(188));      % 1 bit E5a Data validity status (0-valid)

        validWord = 1;

    case 2   %--- Ephemeris (1/3) and GST -------------------------------------
        eph.IODnav2   = bin2dec(navWord(7:16));                           % 10 bit
        eph.M_0       = twosComp2dec(navWord(17:48))   * 2^(-31) * galPi; % 32 bit in rad
        eph.OmegaDot  = twosComp2dec(navWord(49:72))   * 2^(-43) * galPi; % 24 bit in rad
        eph.e         = bin2dec(navWord(73:104))        * 2^(-33);         % 32 bit
        eph.sqrtA     = bin2dec(navWord(105:136))       * 2^(-19);         % 32 bit in m^0.2
        eph.Omega_0   = twosComp2dec(navWord(137:168))   * 2^(-31) * galPi; % 32 bit in rad
        eph.iDot      = twosComp2dec(navWord(169:182)) * 2^(-43) * galPi; % 14 bit in rad

        validWord = 1;

    case 3   %--- Ephemeris (2/3) and GST -------------------------------
        eph.IODnav3   = bin2dec(navWord(7:16));                           % 10 bit
        eph.i_0       = twosComp2dec(navWord(17:48))   * 2^(-31) * galPi; % 32 bit in rad
        eph.omega     = twosComp2dec(navWord(49:80))  * 2^(-31) * galPi; % 32 bit in rad
        eph.deltan    = twosComp2dec(navWord(81:96))   * 2^(-43) * galPi; % 16 bit in rad
        eph.CUC       = twosComp2dec(navWord(97:112))   * 2^(-29);         % 16 bit in rad
        eph.CUS       = twosComp2dec(navWord(113:128))   * 2^(-29);         % 16 bit in rad
        eph.CRC       = twosComp2dec(navWord(129:144))  * 2^(-5);          % 16 bit in m
        eph.CRS       = twosComp2dec(navWord(145:160)) * 2^(-5);          % 16 bit in m
        eph.t_oe      = bin2dec(navWord(161:174))        * 60;              % 14 bit in s

        validWord = 1;

    case 4   %--- Ephemeris (3/3), GST-UTC conversion, GST-GPS conversion and TO W
        eph.IODnav4   = bin2dec(navWord(7:16));                    % 10 bit
        eph.CIC       = twosComp2dec(navWord(17:32))   * 2^(-29);  % 16 bit in rad
        eph.CIS       = twosComp2dec(navWord(33:48))   * 2^(-29);  % 16 bit in rad
        eph.A0        = twosComp2dec(navWord(49:80))    * 2^(-30);  % 32 bit
        eph.A1        = twosComp2dec(navWord(81:104))   * 2^(-50);  % 24 bit
        eph.delt_LS   = twosComp2dec(navWord(105:112));              % 8 bit
        eph.t_ot      = bin2dec(navWord(113:120))        * 3600;     % 8 bit
        eph.WN_ot     = bin2dec(navWord(121:128));                   % 8 bit
        eph.WN_LSF    = bin2dec(navWord(129:136));                   % 8 bit
        eph.DN        = bin2dec(navWord(137:139));                   % 3 bit
        eph.delt_LSF  = twosComp2dec(navWord(140:147));             % 8 bit
        eph.t_og      = bin2dec(navWord(148:155))      * 3600;     % 8 bit
        eph.A0_G      = twosComp2dec(navWord(156:171))  * 2^(-35);  % 16 bit
        eph.A1_G      = twosComp2dec(navWord(172:183)) * 2^(-51);  % 12 bit
        eph.WN_og     = bin2dec(navWord(184:189));                 % 6 bit

        validWord = 1;

        %case 5   %--- Almanac -------------------------------------

        %case 6   %--- Almanac -------------------------------------


end % switch word type

end % all pages


