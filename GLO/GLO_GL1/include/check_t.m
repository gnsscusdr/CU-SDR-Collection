function corrTime = check_t(time)
%CHECK_T accounting for beginning or end of day crossover.
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
%GLONASS modification by Jakob Almqvist
%
% CVS record:
% $Id: check_t.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
%==========================================================================

half_day = 43200;     % seconds

corrTime = time;

if time > half_day
    corrTime = time - 2*half_day;
elseif time < -half_day
    corrTime = time + 2*half_day;
end
%%%%%%% end check_t.m  %%%%%%%%%%%%%%%%%