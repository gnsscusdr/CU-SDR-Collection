function [eph, SOW] = ephemeris(bits, PRN)
%Function decodes ephemerides and SOW from the given bit stream. The stream
%(array) in the parameter BITS must contain 1500 bits. The first element in
%the array must be the first bit of a subframe. The subframe ID of the
%first subframe in the array is not important.
%
%Function does not check parity!
%
%[eph, SOW] = ephemeris(bits, PRN)
%
%   Inputs:
%       bits        - bits of the navigation messages.
%                   Type is character array and it must contain only
%                   characters '0' or '1'.
%       PRN         - PRN number to seperate the decoding process.
%                   GEO & MEO/IGSO have different message structure.
%
%   Outputs:
%       SOW         - Second Of Week (SOW) of the first sub-frame in the bit
%                   stream (in seconds)
%       eph         - SV ephemeris
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Written by Daehee Won, Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

% Modified by Daehee Won.
% Final revision: Sep. 20, 2012.


%--- Decode Navigation mesaage ----------------------
% The task is to select the necessary bits and convert them to decimal
% numbers. For more details on sub-frame contents please refer to
% Beidou ICD (Version 1.0, December 2012).


eph = eph_structure_init();
SOW = inf;
%% Check if there is enough data ==========================================
if length(bits) < 1500
    error('The parameter BITS must contain 1500 bits!');
end

%% Check if the parameters are strings ====================================

% Pi used in the GPS coordinate system
BeidouPi = 3.1415926535898;
% Preamble for Beidou
preamble = [1 1 1 0 0 0 1 0 0 1 0];
Ipreamble = [0 0 0 1 1 1 0 1 1 0 1];

% 'bits' should be row vector for 'bin2dec' function.
[a, b] = size(bits);
if a > b
    bits = bits';
end

if isequal(bits(1:11),preamble) || isequal(bits(1:11),Ipreamble)
    %--- Correct polarity of the data bits in all 10 words ----------------
    if isequal(bits(1:11),Ipreamble)
        bits = bits*-1 +1;
    end
else
    disp(['Preamble does NOT match in for PRN ',num2str(PRN),'!']);
    return;
end

%% Ephemeris for GEO ======================================================
if ((1 <= PRN) && (PRN <=5 )) || ((59 <= PRN) && (PRN <= 63))
    
    % GEO has 10 words in subframe 1. Each word has 150 navigation bits
    % i.e., 300 raw bits (300 ms)
    
    % BCH decoding & interleaving should be apllied!!!
    
    a1_msb       = []; a1_lsb =[];
    C_uc_msb     = []; C_uc_lsb = [];
    e_msb        = []; e_lsb = [];
    C_ic_msb     = []; C_ic_lsb = [];
    i_0_msb      = []; i_0_lsb = [];
    omegaDot_msb = []; omegaDot_lsb = [];
    omega_msb    = []; omega_lsb = [];
    
    D1_valid = inf(1,10);
    
    for i = 1:50
        
        %--- "Cut" one sub-frame's 300 bits -------------------------------
        subframe = bits(300*(i-1)+1 : 300*i);
        
        % Deinterleaving---------------------------------------------------
        % The first word
        DeSubframe = subframe(1:30);
        % Other 4 words
        for k = 1:4 % 150 LSBs are reserved
            DeSubframe = [DeSubframe, ...
                subframe(k*30+1:2:k*30+22), ...
                subframe(k*30+2:2:k*30+22), ...
                subframe(k*30+23:2:k*30+30), ...
                subframe(k*30+24:2:k*30+30)]; %#ok<AGROW>
        end
        subframe = DeSubframe;
        
        %--- Decode the sub-frame id ------------------------------------------
        % Decode the received signal in code using an [15,11] BCH decoder
        [decoded,cnumerr] = bchdec( gf(subframe(16:30),1),15,11 ); 
        if cnumerr ~= -1
            subframe(16:26) = double(decoded.x);
            IDbin = num2str(subframe(16:18));
            subframeID = bin2dec(IDbin);
        else
            disp(['ID of subframe ',num2str(i),'for PRN# ',...
                num2str(PRN),'can not be decoded!']);
            break
        end
        
        % For more details on sub-frame contents please refer to Beidou IS.
        % Only Subframe 1 includes the navigation message.
        if subframeID == 1
            % BCH decode
            for ind = 1:4
                codeword = [subframe(30*ind+1:30*ind+11), subframe(30*ind+23:30*ind+26);...
                    subframe(30*ind+12:30*ind+22), subframe(30*ind+27:30*ind+30)];
                [decoded,cnumerr] = bchdec( gf(codeword,1),15,11 );
                decodeBCH = double(decoded.x);
                if ~isequal(cnumerr,[-1 -1]')
                    subframe(30*ind+1:30*ind+11) = decodeBCH(1,1:11);
                    subframe(30*ind+12:30*ind+22) = decodeBCH(2,1:11);
                else
                    disp(['BCH decoding for PRN#',num2str(PRN),'fails!']);
                    break
                end
            end
            
            subframe = num2str(subframe')';
            Pnum1 = bin2dec(subframe(43:46));   % range: 1-10
            
            if SOW == inf
                SOW = bin2dec(subframe([19:26, 31:42])) - 0.6*(i-1);
                eph.SOW = SOW;
            end
            
            switch Pnum1
                case 1
                    eph.SatH1  = bin2dec( subframe(47) );       % Autonomous Satellite Health flag: 0, 1
                    eph.IODC   = bin2dec( subframe(48:52) );    % Issue of Data, Clock
                    eph.URAI   = bin2dec( subframe(61:64) );    % User Range Accuracy Index: 0~15
                    eph.WN     = bin2dec( subframe(65:77) );    % Week Number: 0~8191
                    
                    eph.t_oc   = bin2dec( subframe([78:82, 91:102]) ) *2^(3);	% Clock Correction Parameters
                    eph.T_GD_1 = twosComp2dec( subframe(103:112) ) * 0.1 * 10^(-9);  % Equipment Group Delay Differential
                    D1_valid(1) = 1;
                    
                case 2  % Ionospheric Delay Model Parameters (alpha, beta)
                    eph.alpha0 = twosComp2dec( subframe([47:52, 61:62]) ) *2^(-30);  % [s]
                    eph.alpha1 = twosComp2dec( subframe(63:70) ) *2^(-27);           % [s/pi]
                    eph.alpha2 = twosComp2dec( subframe(71:78) ) *2^(-24);           % [s/pi^2]
                    eph.alpha3 = twosComp2dec( subframe([79:82, 91:94]) ) *2^(-24);  % [s/pi^3]
                    
                    eph.beta0  = twosComp2dec( subframe(95:102) ) *2^(11);           % [s]
                    eph.beta1  = twosComp2dec( subframe(103:110) ) *2^(14);          % [s/pi]
                    eph.beta2  = twosComp2dec( subframe([111:112, 121:126]) ) *2^(16); % [s/pi^2]
                    eph.beta3  = twosComp2dec( subframe(127:134) ) *2^(16);          % [s/pi^3]
                    
                    D1_valid(2) = 1;
                    
                case 3  % Clock Correction Parameters
                    eph.a0     = twosComp2dec( subframe([101:112, 121:132]) ) *2^(-33);	% [s]
                    a1_msb     = subframe(133:136);
                    D1_valid(3) = 1;
                    
                case 4
                    % Clock Correction Parameters (Cont.)
                    a1_lsb     = subframe([47:52, 61:72]);
                    eph.a2     = twosComp2dec( subframe([73:82, 91]) ) *2^(-66);	% [s/s^2]
                    
                    % Issue of Data, Ephemeris (IODE)
                    eph.IODE   = bin2dec( subframe(92:96) );
                    
                    % Ephemeris Parameters
                    eph.deltan = twosComp2dec( subframe(97:112) ) * 2^(-43) * BeidouPi;      % [pi/s]
                    C_uc_msb   = subframe(121:134);
                    D1_valid(4) = 1;
                    
                case 5  % Ephemeris Parameters (Cont.)
                    C_uc_lsb   = subframe(47:50);
                    eph.M_0    = twosComp2dec( subframe([51:52, 61:82, 91:98]) ) * 2^(-31) * BeidouPi; % [pi]
                    eph.C_us   = twosComp2dec( subframe([99:112, 121:124]) ) * 2^(-31);
                    e_msb      = subframe(125:134);
                    D1_valid(5) = 1;
                    
                case 6  % Ephemeris Parameters (Cont.)
                    e_lsb      = subframe([47:52, 61:76]);
                    eph.sqrtA  = bin2dec( subframe([77:82, 91:112, 121:124]) ) * 2^(-19);
                    C_ic_msb   = subframe(125:134);
                    D1_valid(6) = 1;
                    
                case 7  % Ephemeris Parameters (Cont.)
                    C_ic_lsb   = subframe([47:52, 61:62]);
                    eph.C_is   = twosComp2dec( subframe(63:80) ) * 2^(-31);
                    eph.t_oe   = bin2dec( subframe([81:82, 91:105]) )  * 2^3;
                    i_0_msb    = subframe([106:112, 121:134]);
                    D1_valid(7) = 1;
                    
                case 8  % Ephemeris Parameters (Cont.)
                    i_0_lsb    = subframe([47:52, 61:65]);
                    eph.C_rc   = twosComp2dec( subframe([66:82, 91]) ) * 2^(-6);
                    eph.C_rs   = twosComp2dec( subframe(92:109) ) * 2^(-6);
                    omegaDot_msb = subframe([110:112, 121:136]);
                    D1_valid(8) = 1;
                    
                case 9  % Ephemeris Parameters (Cont.)
                    omegaDot_lsb = subframe(47:51);
                    eph.omega_0  = twosComp2dec( subframe([52, 61:82, 91:99]) ) * 2^(-31) * BeidouPi;
                    omega_msb    = subframe([100:112, 121:134]);
                    D1_valid(9) = 1;
                case 10  % Ephemeris Parameters (Cont.)
                    omega_lsb    = subframe(47:51);
                    eph.iDot     = twosComp2dec( subframe([52, 61:73]) ) * 2^(-43) * BeidouPi;
                    D1_valid(10) = 1;
            end
        end % if subframeID == 1
        
        %% MSB & LSB combination
        if length([a1_msb a1_lsb]) == 22
            eph.a1     = twosComp2dec( [a1_msb, a1_lsb] ) *2^(-50);
        end
        
        if length([C_uc_msb, C_uc_lsb]) == 18
            eph.C_uc   = twosComp2dec( [C_uc_msb, C_uc_lsb] ) * 2^(-31);
        end
        
        if length([e_msb, e_lsb]) == 32
            eph.e      = bin2dec( [e_msb, e_lsb] ) * 2^(-33);
        end
        
        if length([C_ic_msb, C_ic_lsb]) == 18
            eph.C_ic   = twosComp2dec( [C_ic_msb, C_ic_lsb] ) * 2^(-31);
        end
        
        if length([i_0_msb, i_0_lsb]) == 32
            eph.i_0    = twosComp2dec( [i_0_msb, i_0_lsb] ) * 2^(-31) * BeidouPi;
        end
        
        if length([omegaDot_msb, omegaDot_lsb]) == 24
            eph.omegaDot = twosComp2dec( [omegaDot_msb, omegaDot_lsb] ) * 2^(-43) * BeidouPi;
        end
        
        if length([omega_msb, omega_lsb]) == 32
            eph.omega    = twosComp2dec( [omega_msb, omega_lsb] ) * 2^(-31) * BeidouPi;
        end
    end % for i = 1:50
    
    if sum(D1_valid) == 10
        eph.flag = 1;
    end
    
    % Compute the second of week (SOW) of the first sub-frames in the array ====
    % Also correct the SOW. The transmitted SOW is actual SOW of the next
    % subframe and we need the SOW of the first subframe in this data block
    % (the variable subframe at this point contains bits of the last subframe).
    % D2 subframe is 3 seconds long.
    %     subframe = num2str(subframe')';
    %     SOW = bin2dec(subframe([19:26, 31:42]))-27;
    
    %% Ephemeris for MEO/IGSO  ================================================
elseif (6 <= PRN) && (PRN <= 58)
    % Decode all 5 sub-frames =============================================
    
    D1_valid = inf(1,3);
    
    % BCH decoding & interleaving should be apllied!!!
    
    t_oe_msb = [];
    t_oe_lsb = [];
    
    for i = 1:5
        
        %--- "Cut" one sub-frame's bits ---------------------------------------
        subframe = bits(300*(i-1)+1 : 300*i);
        
        % Deinterleaving --------------------------------------------------
        DeSubframe = subframe(1:30);
        for k = 1:9
            DeSubframe = [DeSubframe, ...
                subframe(k*30+1:2:k*30+22), ...
                subframe(k*30+2:2:k*30+22), ...
                subframe(k*30+23:2:k*30+30), ...
                subframe(k*30+24:2:k*30+30)];
        end
        subframe = DeSubframe;
        
        %--- Decode the sub-frame id ------------------------------------------
        % For more details on sub-frame contents please refer to GPS IS.
        %--- Decode the sub-frame id ------------------------------------------
        [decoded,cnumerr] = bchdec( gf(subframe(16:30),1),15,11 );
        if cnumerr ~= -1
            subframe(16:26) = double(decoded.x);
            IDbin = num2str(subframe(16:18));
            subframeID = bin2dec(IDbin);
        else
            disp(['ID of subframe ',num2str(i),'for PRN# ',...
                num2str(PRN),'can not be decoded!']);
            break
        end
        %         subframeID = bin2dec(subframe(16:18));
        
        % do BCH decoding
        if subframeID == 1 || subframeID == 2 || subframeID == 3
            for ind = 1:4
                codeword = [subframe(30*ind+1:30*ind+11), subframe(30*ind+23:30*ind+26);...
                    subframe(30*ind+12:30*ind+22), subframe(30*ind+27:30*ind+30)];
                [decoded,cnumerr] = bchdec( gf(codeword,1),15,11 );
                decodeBCH = double(decoded.x);
                if ~isequal(cnumerr,[-1 -1]')
                    subframe(30*ind+1:30*ind+11) = decodeBCH(1,1:11);
                    subframe(30*ind+12:30*ind+22) = decodeBCH(2,1:11);
                else
                    disp(['BCH decoding for PRN#',num2str(PRN),'fails!']);
                    break
                end
            end
            
            subframe = num2str(subframe')';
            
            if SOW==inf
                SOW = bin2dec(subframe([19:26, 31:42])) - (i-1)*6;
                eph.SOW = SOW;
            end
        end
        
        switch subframeID
            case 1  %--- It is subframe 1 -------------------------------------
                
                eph.SatH1  = bin2dec( subframe(43) );       % Autonomous Satellite Health flag: 0, 1
                eph.IODC   = bin2dec( subframe(44:48) );    % Issue of Data, Clock
                eph.URAI   = bin2dec( subframe(49:52) );    % User Range Accuracy Index: 0~15
                eph.WN     = bin2dec( subframe(61:73) );    % Week Number: 0~8191
                
                eph.t_oc   = bin2dec( subframe([74:82, 91:98]) ) * 2^(3);	% Clock Correction Parameters
                eph.T_GD_1 = twosComp2dec( subframe(99:108) ) * 0.1 * 10^(-9);% Equipment Group Delay Differential
                
                % Ionospheric Delay Model Parameters (alpha, beta)
                eph.alpha0 = twosComp2dec( subframe(127:134) ) *2^(-30); % [s]
                eph.alpha1 = twosComp2dec( subframe(135:142) ) *2^(-27); % [s/pi]
                eph.alpha2 = twosComp2dec( subframe(151:158) ) *2^(-24); % [s/pi^2]
                eph.alpha3 = twosComp2dec( subframe(159:166) ) *2^(-24); % [s/pi^3]
                
                eph.beta0  = twosComp2dec( subframe([167:172,181:182]) ) *2^(11); % [s]
                eph.beta1  = twosComp2dec( subframe(183:190) ) *2^(14);           % [s/pi]
                eph.beta2  = twosComp2dec( subframe(191:198) ) *2^(16);           % [s/pi^2]
                eph.beta3  = twosComp2dec( subframe([199:202, 211:214]) ) *2^(16);% [s/pi^3]
                
                % Clock Correction Parameters
                eph.a2     = twosComp2dec( subframe(215:225) ) *2^(-66);            % [s/s^2]
                eph.a0     = twosComp2dec( subframe([226:232, 241:257]) ) *2^(-33);	% [s]
                eph.a1     = twosComp2dec( subframe([258:262, 271:287]) ) *2^(-50);	% [s/s]
                
                % Issue of Data, Ephemeris (IODE)
                eph.IODE   = bin2dec( subframe(288:292) );
                
                D1_valid(1) = 1;
                
                
            case 2  %--- It is subframe 2 -------------------------------------
                % Ephemeris Parameters
                eph.deltan = twosComp2dec( subframe([43:52, 61:66])) * 2^(-43) * BeidouPi;
                eph.C_uc   = twosComp2dec( subframe([67:82,  91:92])) * 2^(-31);
                eph.M_0    = twosComp2dec( subframe([93:112, 121:132])) * 2^(-31) * BeidouPi;
                eph.e      = bin2dec( subframe([133:142, 151:172])) * 2^(-33);
                
                eph.C_us   = twosComp2dec( subframe(181:198)) * 2^(-31);
                eph.C_rc   = twosComp2dec( subframe([199:202, 211:224])) * 2^(-6);
                eph.C_rs   = twosComp2dec( subframe([225:232, 241:250])) * 2^(-6);
                
                eph.sqrtA  = bin2dec( subframe([251:262, 271:290])) * 2^(-19);
                t_oe_msb   = subframe(291:292);
                
                D1_valid(2) = 1;
                
            case 3  %--- It is subframe 3 -------------------------------------
                t_oe_lsb     = subframe([43:52, 61:65]);
                eph.i_0      = twosComp2dec( subframe([66:82,    91:105])) * 2^(-31) * BeidouPi;
                eph.C_ic     = twosComp2dec( subframe([106:112, 121:131])) * 2^(-31);
                eph.omegaDot = twosComp2dec( subframe([132:142, 151:163])) * 2^(-43) * BeidouPi;
                eph.C_is     = twosComp2dec( subframe([164:172, 181:189])) * 2^(-31);
                eph.iDot     = twosComp2dec( subframe([190:202, 211]    )) * 2^(-43) * BeidouPi;
                eph.omega_0  = twosComp2dec( subframe([212:232, 241:251])) * 2^(-31) * BeidouPi;
                eph.omega    = twosComp2dec( subframe([252:262, 271:291])) * 2^(-31) * BeidouPi;
                
                D1_valid(3) = 1;
                
            case 4  %--- It is subframe 4 -------------------------------------
                % Not decoded at the moment.
                
            case 5  %--- It is subframe 5 -------------------------------------
                % Not decoded at the moment.
                
        end % switch subframeID ...
        
    end % for all 5 sub-frames ...
    
    % Ephemeris reference time
    if length([t_oe_msb, t_oe_lsb]) == 17
        eph.t_oe = bin2dec([t_oe_msb, t_oe_lsb]) * 2^3;
    else
        eph.t_oe = nan;
    end
    
    if sum(D1_valid) == 3
        eph.flag = 1;
    end
    % Compute the second of week (SOW) of the first sub-frames in the array ====
    % Also correct the SOW. The transmitted SOW is actual SOW of the next
    % subframe and we need the SOW of the first subframe in this data block
    % (the variable subframe at this point contains bits of the last subframe).
    % D1 subframe is 6 seconds long.
else
    disp('PRN is NOT in range between 1-63 ');
    
end
