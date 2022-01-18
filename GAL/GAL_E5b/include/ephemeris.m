function [eph] = ephemeris(navBits,eph)
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
%--------------------------------------------------------------------------

%% Check if there is enough data ==========================================
% This is 15 full pages (one subframe)
if length(navBits) < 7500
    error('The ephemeris bitstream must contain 7500 symbols!');
end

%--- Check polarity of the data bits ----------------------------------
sync_bits = [0 1 0 1 1 0 0 0 0 0];

% Convert CNAV-producing convolutional code polynomials to trellis description
% Note that the difference from GPS is that the second branch G2 is
% inverted at the end (see ICD)
trellis = poly2trellis(7,[171 ~133]);

%Creates a cyclic redundancy code (CRC) detector System object
crcDet = comm.CRCDetector([24 23 18 17 14 11 10 7 6 5 4 3 1 0]);

% Viterbi traceback depth for vitdec(function)
tblen = 35;

%Indicator the requisite messages are all decoded
word_valid = inf(1,6);
%% Define required constants ==============================================
% Pi used in the Galileo coordinate system (same as for GPS)
galPi = 3.1415926535898;
eph = eph_structure_init();

%% Decode Messages ========================================================
for ii = 1:2:30
    
    %--- Pull out bits for both page parts --------------------------------
    pagePart1 = navBits(250*(ii-1)+1 : 250* ii   )';
    pagePart2 = navBits(250* ii   +1 : 250*(ii+1))';
    
    %--- Correct polarity of the all data bits according to preamble bits
    if(~isequal(pagePart1(1:10),sync_bits))
        pageSymInt1 = not(pagePart1(11:250));
        pageSymInt2 = not(pagePart2(11:250));
    else
        pageSymInt1 = pagePart1(11:250);
        pageSymInt2 = pagePart2(11:250);
    end
    
    %--- De-interleave page parts -----------------------------------------
    symMat1 = reshape(pageSymInt1,30,8)';
    pageSym1 = reshape(symMat1,1,[])';
    symMat2 = reshape(pageSymInt2,30,8)';
    pageSym2 = reshape(symMat2,1,[])';
    
    
    %--- Remove convolutional encoding from page parts --------------------
    decBits1 = vitdec(pageSym1,trellis,tblen,'trunc','hard');
    decBits2 = vitdec(pageSym2,trellis,tblen,'trunc','hard');
    
    %--- Reconstruct full page --------------------------------------------
    if (decBits1(1) == 0) && (decBits2(1) == 1)
        page = [decBits1(1:114)' decBits2(1:106)'];
        part = 1;
    elseif (decBits1(1) == 1) && (decBits2(1) == 0)
        page = [decBits2(1:114)' decBits1(1:106)'];
        part = 2;
    else
        break % Pages are not even and odd, assume sync pattern is off
    end
    
    %--- Check the CRC ----------------------------------------------------
    [~,frmError] = step(crcDet,page');
    if (frmError)
        continue
    end
    
    %--- Pull out navigation data word ------------------------------------
    navWordDec = [page(3:114) page(117:132)];
    
    %--- Convert navigation data word to binary ---------------------------
    navWord = dec2bin(navWordDec)';
    
    %--- Decode the message type ------------------------------------------
    wordType = bin2dec(navWord(1:6));
    
    
    %--- Decode sub-frame based on the message type -----------------------
    switch wordType
        case 1   %--- Ephemeris (1/4) -------------------------------------
            eph.IODnav1   = bin2dec(navWord(7:16));
            eph.t_oe      = bin2dec(navWord(17:30))        * 60;    % in s
            eph.M_0       = twosComp2dec(navWord(31:62))   * 2^(-31) * galPi; % in rad
            eph.e         = bin2dec(navWord(63:94))        * 2^(-33);
            eph.sqrtA     = bin2dec(navWord(95:126))       * 2^(-19);  % in m^0.2
            word_valid(1) = 1;
            
        case 2   %--- Ephemeris (2/4) -------------------------------------
            eph.IODnav2   = bin2dec(navWord(7:16));
            eph.Omega_0   = twosComp2dec(navWord(17:48))   * 2^(-31) * galPi; % in rad
            eph.i_0       = twosComp2dec(navWord(49:80))   * 2^(-31) * galPi; % in rad
            eph.omega     = twosComp2dec(navWord(81:112))  * 2^(-31) * galPi; % in rad
            eph.iDot      = twosComp2dec(navWord(113:126)) * 2^(-43) * galPi; % in rad
            word_valid(2) = 1;
            
        case 3   %--- Ephemeris (3/4), SISA -------------------------------
            eph.IODnav3   = bin2dec(navWord(7:16));
            eph.OmegaDot  = twosComp2dec(navWord(17:40))   * 2^(-43) * galPi; % in rad
            eph.deltan    = twosComp2dec(navWord(41:56))   * 2^(-43) * galPi; % in rad
            eph.CUC       = twosComp2dec(navWord(57:72))   * 2^(-29);         % in rad
            eph.CUS       = twosComp2dec(navWord(73:88))   * 2^(-29);         % in rad
            eph.CRC       = twosComp2dec(navWord(89:104))  * 2^(-5);          % in m
            eph.CRS       = twosComp2dec(navWord(105:120)) * 2^(-5);          % in m
            %eph.SISA      = navWord(121:128) % content not currently defined
            word_valid(3) = 1;
            
        case 4   %--- SVID, Ephemeris (4/4), Clock Correction -------------
            eph.IODnav4   = bin2dec(navWord(7:16));
            eph.SVID      = bin2dec(navWord(17:22));
            eph.CIC       = twosComp2dec(navWord(23:38))   * 2^(-29); % in rad
            eph.CIS       = twosComp2dec(navWord(39:54))   * 2^(-29); % in rad
            eph.t_oc      = bin2dec(navWord(55:68))        * 60;      % in s
            eph.a_f0      = twosComp2dec(navWord(69:99))   * 2^(-34); % in s
            eph.a_f1      = twosComp2dec(navWord(100:120)) * 2^(-46); % in s/s
            eph.a_f2      = twosComp2dec(navWord(121:126)) * 2^(-59); % in s/(s^2)
            %             eph.a_f2      = 0; % --------------------------------------------------------------------------
            word_valid(4) = 1;
            
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
            eph.TOW = eph.TOW - ii + 1;
            if part == 2
                eph.TOW = eph.TOW + 1;
            end
            word_valid(5) = 1;
            
        case 6   %--- GST-UTC Conversion ----------------------------------
            eph.A0        = twosComp2dec(navWord(7:38))    * 2^(-30);
            eph.A1        = twosComp2dec(navWord(39:62))   * 2^(-50);
            eph.delt_LS   = twosComp2dec(navWord(63:70));
            eph.t_ot      = bin2dec(navWord(71:78))        * 3600;
            eph.WN_ot     = bin2dec(navWord(79:86));
            eph.WN_LSF    = bin2dec(navWord(87:94));
            eph.DN        = bin2dec(navWord(95:97));
            eph.delt_LSF  = twosComp2dec(navWord(98:105));
            word_valid(6) = 1;
            
            %case 7   %--- Almanac -------------------------------------
            
            %case 8   %--- Almanac -------------------------------------
            
            %case 9   %--- Almanac -------------------------------------
            
        case 10  %--- GPS-GST Conversion ---------------------------------
            eph.A0_G      = twosComp2dec(navWord(87:102))  * 2^(-35);
            eph.A1_G      = twosComp2dec(navWord(103:114)) * 2^(-51);
            eph.t_og      = bin2dec(navWord(115:122))      * 3600;
            eph.WN_og     = bin2dec(navWord(123:128));
            
    end % switch word type
    if sum(word_valid) == 6
        % Check if the Issue of Data (IOD ) values are the same
        IODnav = [eph.IODnav1 eph.IODnav2 eph.IODnav3 eph.IODnav4];
        if (length(unique(IODnav)) == 1)
            eph.flag = 1;
        else
            disp('    The IOD values of each pages are different!')
        end
        
        break
    end
    
end % all pages
