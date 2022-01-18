function [status, C] = dataVerification(data)
% This function uses the data verification algorithm described in the
% GLONASS ICD. It does not correct the data, it only checks if it is
% correct or not.
%
%status = dataVerification(data)
%
%   Inputs: 
%       data        - 85 bits. 77 Navigation bits and 8 Hamming
%                   code checking bits.
%
%   Outputs: 
%       status      - returns 1 if correct or 0 if not correct.
%

%--------------------------------------------------------------------------
%
%                       Written by: Jakob Almqvist
%
%--------------------------------------------------------------------------
index1 = [9 10 12 13 15 17 19 20 22 24 26 28 30 32 34 35 37 39 41 ...
 43 45 47 49 51 53 55 57 59 61 63 65 66 68 70 72 74 76 78 80 82 84];
C(1,1) = xor(data(1,1),mod(sum(data(1,index1)),2));

index2 = [9 11 12 14 15 18 19 21 22 25 26 29 30 33 34 36 37 40 41 ...
 44 45 48 49 52 53 56 57 60 61 64 65 67 68 71 72 75 76 79 80 83 84];
C(1,2) = xor(data(1,2),mod(sum(data(1,index2)),2));

index3 = [10:12 16:19 23:26 31:34 38:41 46:49 54:57 62:65 69:72 77:80 85];
C(1,3) = xor(data(1,3),mod(sum(data(1,index3)),2));

index4 = [13:19 27:34 42:49 58:65 73:80];
C(1,4) = xor(data(1,4),mod(sum(data(1,index4)),2));

index5 = [20:34 50:65 81:85];
C(1,5) = xor(data(1,5),mod(sum(data(1,index5)),2));

index6 = 35:65;
C(1,6) = xor(data(1,6),mod(sum(data(1,index6)),2));

C(1,7) = xor(data(1,7),mod(sum(data(1,66:85)),2));

C(1,8) = xor(mod(sum(data(1:8)),2),mod(sum(data(9:85)),2));

 if isempty(find(C,1)) || (length(find(C(1,1:7))) == 1 && C(1,8) == 1)
     status = 1;
 else
     status = 0;
 end