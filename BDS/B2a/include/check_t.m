function corrTime = check_t(time)
%CHECK_T accounting for beginning or end of week crossover.
%
%corrTime = check_t(time);
%
%   Inputs:
%       time        - time in seconds
%
%   Outputs:
%       corrTime    - corrected time (seconds)

%Kai Borre 04-01-96
%Copyright (c) by Kai Borre
%
% CVS record:
% $Id: check_t.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
%==========================================================================

half_week = 302400;     % seconds

corrTime = time;
% The message type 10, 11, and 30-37 data sets are transmitted by the SV for 
% periods of two hours.So the following condition indicates that beginning 
% or end of week crossover has occurred, so the time should be corrected!
% Refer to IS-GPS-200H P163.
if time > half_week
    corrTime = time - 2*half_week;
elseif time < -half_week
    corrTime = time + 2*half_week;
end
%%%%%%% end check_t.m  %%%%%%%%%%%%%%%%%