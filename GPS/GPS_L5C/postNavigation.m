function [navSolutions, eph] = postNavigation(trackResults, settings)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original work by Darius Plausinaitis,Peter Rinder, 
% Nicolaj Bertelsen and Dennis M. Akos
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
%$Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $

%% Check is there enough data to obtain any navigation solution ===========
% It is necessary to have at least three messages (type 10, 11 and 
% anyone of 30-37) to find satellite coordinates. Then receiver 
% position can be found. The function requires at least 3 message.
% One message length is 6 seconds, therefore we need at least 18 sec long
% record (3 * 6 = 18 sec = 18000ms). We add extra seconds for the cases,
% when tracking has started in a middle of a message.Therefor, at least
% 24000ms signal is required to decode requisite ephemeris for fix.
if (settings.msToProcess < 24000)
    % Show the error message and exit
    disp('Record is to short. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Pre-allocate space =======================================================
% Starting positions of the first message in the input bit stream 
% trackResults.I_P in each channel. The position is PRN code count
% since start of tracking. Corresponding value will be set to inf 
% if no valid preambles were detected in the channel.
subFrameStart  = inf(1, settings.numberOfChannels);

% Time Of Week (TOW) of the first message(in seconds). Corresponding value
% will be set to inf if no valid preambles were detected in the channel.
TOW  = inf(1, settings.numberOfChannels);

%--- Make a list of channels excluding not tracking channels ---------------
activeChnList = find([trackResults.status] ~= '-');

%% Decode ephemerides =======================================================
for channelNr = activeChnList
    
    % Get PRN of current channel
    PRN = trackResults(channelNr).PRN;
    
    fprintf('Decoding CNAV for PRN %02d -------------------- \n', PRN);

    %=== Decode ephemerides and TOW of the first sub-frame ==================
    [eph(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
                                  CNAVdecoding(trackResults(channelNr).I_P);  %#ok<AGROW>

    %--- Exclude satellite if it does not have the necessary cnav data ------
    if (eph(PRN).idValid(1) ~= 10 || eph(PRN).idValid(2) ~= 11 ...
        || ~sum(eph(PRN).idValid(3:10) == (30:37)) )

        %--- Exclude channel from the list ----------------------------------
        activeChnList = setdiff(activeChnList, channelNr);
        
        %--- Print CNAV decoding information for current PRN ----------------
        if (eph(PRN).idValid(1) ~= 10)
            fprintf('    Message type 10 for PRN %02d not decoded.\n', PRN);
        end
        
        if (eph(PRN).idValid(2) ~= 11)
            fprintf('    Message type 11 for PRN %02d not decoded.\n', PRN);
        end
        
        if (~sum(eph(PRN).idValid(3:10) == (30:37)))
            fprintf('    None of message type 30-37 for PRN %02d decoded.\n', PRN);
        end
        fprintf('    Channel for PRN %02d excluded!!\n', PRN);
    
    else
        fprintf('    Three requisite messages for PRN %02d all decoded!\n', PRN);
    end 
end %  channelNr = activeChnList

%% Check if the number of satellites is still above 3 =====================
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Set measurement-time point and step  =====================================
% Find start and end of measurement point locations in IF signal stream with available
% measurements
sampleStart = zeros(1, settings.numberOfChannels);
sampleEnd = inf(1, settings.numberOfChannels);
for channelNr = activeChnList
    sampleStart(channelNr) = ...
        trackResults(channelNr).absoluteSample(subFrameStart(channelNr));
    sampleEnd(channelNr) = trackResults(channelNr).absoluteSample(end);
end

% Second term is to make space to aviod index exceeds matrix dimensions, 
% thus a margin of 1 is added.
sampleStart = max(sampleStart) + 1;
sampleEnd = min(sampleEnd) - 1;
 
%--- Measurement step in unit of IF samples -------------------------------
measSampleStep = fix(settings.samplingFreq * settings.navSolPeriod/1000);

%---  Number of measurment point from measurment start to end ------------- 
measNrSum = fix((sampleEnd-sampleStart)/measSampleStep);

%% Initialization =========================================================
% Set the satellite elevations array to INF to include all satellites for
% the first calculation of receiver position. There is no reference point
% to find the elevation angle as there is no receiver position estimate at
% this point.
satElev  = inf(1, settings.numberOfChannels);

% Save the active channel list. The list contains satellites that are
% tracked and have the required ephemeris data. In the next step the list
% will depend on each satellite's elevation angle, which will change over
% time.  
readyChnList = activeChnList;

% Set local time to inf for first calculation of receiver position. After
% first fix, localTime will be updated by measurement sample step.
localTime = inf;

%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################

fprintf('Positions are being computed. Please wait... \n');
for currMeasNr = 1:measNrSum
   
    fprintf('Fix: Processing %02d of %02d \n', currMeasNr,measNrSum);
    
    %% Initialization of current measurement ==============================          
    % Exclude satellites, that are belove elevation mask 
    activeChnList = intersect(find(satElev >= settings.elevationMask), ...
                              readyChnList);

    % Save list of satellites used for position calculation
    navSolutions.PRN(activeChnList, currMeasNr) = ...
                                        [trackResults(activeChnList).PRN]; 

    % These two lines help the skyPlot function. The satellites excluded
    % do to elevation mask will not "jump" to possition (0,0) in the sky
    % plot.
    navSolutions.el(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions.az(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
                                     
    % Signal transmitting time of each channel at measurement sample location
    navSolutions.transmitTime(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);
    navSolutions.satClkCorr(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);                                                                  
    % Position index of current measurement time in IF signal stream
    % (in unit IF signal sample point)
    currMeasSample = sampleStart + measSampleStep*(currMeasNr-1);
                                                                      
%% Find pseudoranges ======================================================
    % Raw pseudorange = (localTime - transmitTime) * light speed (in m)
    % All output are 1 by settings.numberOfChannels columme vecters.
    [navSolutions.rawP(:, currMeasNr),transmitTime,localTime]=  ...
                     calculatePseudoranges(trackResults,subFrameStart,TOW, ...
                     currMeasSample,localTime,activeChnList, settings);     
    % Save transmitTime
    navSolutions.transmitTime(activeChnList, currMeasNr) = ...
                                        transmitTime(activeChnList);

%% Find satellites positions and clocks corrections =======================
    % Outputs are all colume vectors corresponding to activeChnList
    [satPositions, satClkCorr] = satpos(transmitTime(activeChnList), ...
                                        [trackResults(activeChnList).PRN], eph); 
                                                                      
    % Save satClkCorr
    navSolutions.satClkCorr(activeChnList, currMeasNr) = satClkCorr;
%% Find receiver position =================================================
    % 3D receiver position can be found only if signals from more than 3
    % satellites are available  
    if size(activeChnList, 2) > 3

        %=== Calculate receiver position ==================================
        % Correct pseudorange for SV clock error
        clkCorrRawP = navSolutions.rawP(activeChnList, currMeasNr)' + ...
                                                   satClkCorr * settings.c;

        % Calculate receiver position
        [xyzdt,navSolutions.el(activeChnList, currMeasNr), ...
         navSolutions.az(activeChnList, currMeasNr), ...
         navSolutions.DOP(:, currMeasNr)] =...
                                 leastSquarePos(satPositions, clkCorrRawP,settings);

        %=== Save results ===========================================================
        % Receiver position in ECEF
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);       
		% For first calculation of solution, clock error will be set 
        % to be zero
        if (currMeasNr == 1)
        navSolutions.dt(currMeasNr) = 0;  % in unit of (m)
        else
            navSolutions.dt(currMeasNr) = xyzdt(4);  
        end
		%=== Correct local time by clock error estimation =================
        localTime = localTime - xyzdt(4)/settings.c;       
        navSolutions.localTime(currMeasNr) = localTime;
        
        % Save current measurement sample location 
        navSolutions.currMeasSample(currMeasNr) = currMeasSample;
        % Update the satellites elevations vector
        satElev = navSolutions.el(:, currMeasNr)';

        %=== Correct pseudorange measurements for clocks errors ===========
        navSolutions.correctedP(activeChnList, currMeasNr) = ...
                navSolutions.rawP(activeChnList, currMeasNr) + ...
                satClkCorr' * settings.c - xyzdt(4);
            
%% Coordinate conversion ==================================================

        %=== Convert to geodetic coordinates ==============================
        [navSolutions.latitude(currMeasNr), ...
         navSolutions.longitude(currMeasNr), ...
         navSolutions.height(currMeasNr)] = cart2geo(...
                                            navSolutions.X(currMeasNr), ...
                                            navSolutions.Y(currMeasNr), ...
                                            navSolutions.Z(currMeasNr), ...
                                            5);
      
        %=== Convert to UTM coordinate system =============================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
                                           navSolutions.longitude(currMeasNr));
        
        % Position in ENU
        [navSolutions.E(currMeasNr), ...
         navSolutions.N(currMeasNr), ...
         navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
                                                xyzdt(3), ...
                                                navSolutions.utmZone);
        
    else
        %--- There are not enough satellites to find 3D position ----------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
                       ': Not enough information for position solution.']);

        %--- Set the missing solutions to NaN. These results will be
        %excluded automatically in all plots. For DOP it is easier to use
        %zeros. NaN values might need to be excluded from results in some
        %of further processing to obtain correct results.
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;

        navSolutions.az(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));
        navSolutions.el(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));

        % TODO: Know issue. Satellite positions are not updated if the
        % satellites are excluded do to elevation mask. Therefore rasing
        % satellites will be not included even if they will be above
        % elevation mask at some point. This would be a good place to
        % update positions of the excluded satellites.
        
        disp('   Exit Program');
        return;


    end % if size(activeChnList, 2) > 3

    %=== Update local time by measurement  step  ====================================
    localTime = localTime + measSampleStep/settings.samplingFreq ;

end %for currMeasNr...
