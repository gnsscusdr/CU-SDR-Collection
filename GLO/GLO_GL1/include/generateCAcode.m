function sampledCode = generateCAcode(PRN,sampFreq,numSamples)

% Function to generate any of the 1-37 (GPS) C/A codes, 120-158 (GPS-SBAS)
% C/A codes, or GLONASS code
%
% Returns 1 for binary 0 and -1 for binary 1
%
% Input
% PRN           - sat id number (GLONASS = 0; PRN = 1-37,120-158)
% sampFreq      - sampling frequency
% numSamples    - number of samples to create
%
% Output
% code          - vector containing the 1 bit/time code sequence (one period)
% sampledCode   - vector containing the sampled code sequence

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (c) Written by Jakob Almqvist,
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


if ( PRN ~= 0 ) % GPS

    codeLength  = 1023;
    codeRate    = 1.023e6;

    if ( PRN > 37 )
        PRN = PRN - 82; % e.g. 120 -> 38
    end

    % the g2s vector holds the appropriate shift of the g2 code to generate
    % the C/A code (ex. for SV#19 - use a G2 shift of g2shift(19,1)=471)
    g2s = [5;6;7;8;17;18;139;140;141;251;252;254;255;256;257;258;469;470;471; ...
        472;473;474;509;512;513;514;515;516;859;860;861;862;863;950;947;948;950; ... % 1-37 (GPS)
        145;175;52;21;237;235;886;657;634;762;355;1012;176; ...
        603;130;359;595;68;386;797;456;499;883;307;127;211; ...
        121;118;163;628;853;484;289;811;202;1021;463;568;904];    % 38-76 [120-158] (GPS-SBAS)

    g2shift=g2s(PRN,1);

    % Generate G1 code
    %   load shift register
    reg = -1*ones(1,10);
    %
    for i = 1:codeLength,
        g1(i) = reg(10);
        save1 = reg(3)*reg(10);     % XOR
        reg(1,2:10) = reg(1:1:9);
        reg(1) = save1;
    end
    %
    % Generate G2 code
    %
    %   load shift register
    reg = -1*ones(1,10);
    %
    for i = 1:codeLength,
        g2(i) = reg(10);    % output from 10th stage
        save2 = reg(2)*reg(3)*reg(6)*reg(8)*reg(9)*reg(10); % XOR
        reg(1,2:10) = reg(1:1:9);
        reg(1) = save2;
    end
    %
    %    Shift G2 code
    %
    g2tmp(1,1:g2shift)=g2(1,codeLength-g2shift+1:1023);
    g2tmp(1,g2shift+1:1023)=g2(1,1:codeLength-g2shift);
    %
    g2 = g2tmp;
    %
    %  Form single sample C/A code by multiplying G1 and G2 point by point
    %
    code = g1.*g2;



else    % GLONASS

    codeLength = 511;
    codeRate = 511e3;

    % Generate code
    %   load shift register
    reg = -1*ones(1,9);
    %
    for i = 1:codeLength,
        code(i) = reg(7);       % output from 7th stage
        save1 = reg(5)*reg(9);  % XOR
        reg(1,2:9) = reg(1:1:8);
        reg(1) = save1;
    end

end

stepSize = codeRate / sampFreq;        % step size (in chips) for sampling

% sample vector at stepSize increments
samples = 0:stepSize:(numSamples*stepSize)-stepSize;

samples = floor(samples);                   % convert to integer values
samples = rem(samples,codeLength);         % handle samples > codeLength
samples = samples + 1;                      % convert from 0-based to 1-based
sampledCode = code(samples);               % grab values from original ranging code
