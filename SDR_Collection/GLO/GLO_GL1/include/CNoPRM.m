function [CNo]= CNoPRM(I,Q,K,M,T)
%Calculate CNo using the Power Ratio Method
%
%[CNo]= CNoVSM(I,Q)
%
%   Inputs:
%       I           - Prompt In Phase values of the signal from Tracking
%       Q           - Prompt Quadrature Phase values of the signal from Tracking
%       K           - Integration time in ms to compute C/No
%       M           - No. of Samples occuring in time T
%       T           - Accumulation Time for Measurement in Tracking
%   Outputs:
%       CNo         - Estimated C/No for the given values of I and Q
%
%--------------------------------------------------------------------------
% Copyright (C) D.M.Akos
% Written by Sirish Jetti
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

%Trim I and Q vectors
I=I(1:20*floor(numel(I)/20));
Q=Q(1:20*floor(numel(Q)/20));
%Reshape I and Q into M columns
I=reshape(I,M,numel(I)/M)';
Q=reshape(Q,M,numel(I)/M)';
%Calculate Wide Band Power
WBP= sum(I.^2+Q.^2,2);

if (WBP ~=0)

    %Calculate Narrow Band Power
    NBP = sum(I,2).^2 + sum(Q,2).^2;
    %Calculate the Noise Power
    NP=NBP./WBP;
    %Calculate the Expected Noise Power
    %The expected noise power value (MuNP ) is the average noise power over h = K/M
    h=K/M;
    %Trim NP vector (Check if this step is required)
    NP=NP(1:h*floor(numel(NP)/h));
    %Reshape NP vector to compute the mean
    NP=reshape(NP,h,numel(NP)/h)';
    %Calculate the expected noise power value (MuNP )
    MuNP=mean(NP,2);
    %Calculate C/No
    CNo=(10*log10(abs((1/T).*(MuNP-1)./(M-MuNP))))';
else
    CNo=zeros(1,floor(numel(I)/K));
end