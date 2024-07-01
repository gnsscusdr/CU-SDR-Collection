function [trackResults, channel]= NB_tracking(fid, channel, settings)
% Performs BDS-3 B1C code and carrier tracking for all channels using
% narrowband correlating approach.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B1C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
%
% Reference: Li, Y., Shivaramaiah, N.C. & Akos, D.M. Design and
% implementation of an open-source BDS-3 B1C/B2a SDR receiver.
% GPS Solut (2019) 23: 60. https://doi.org/10.1007/s10291-019-0853-z
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
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% Number of tracking loop updating
NumToProcess =  round(settings.msToProcess/1000/settings.intTime);

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, NumToProcess);

% Freq of the PRN code:
trackResults.codeFreq       = inf(1, NumToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, NumToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, NumToProcess);
trackResults.I_E            = zeros(1, NumToProcess);
trackResults.I_L            = zeros(1, NumToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, NumToProcess);
trackResults.Q_P            = zeros(1, NumToProcess);
trackResults.Q_L            = zeros(1, NumToProcess);

% for pilot signal
trackResults.Pilot_I_P  = zeros(1, NumToProcess);
trackResults.Pilot_Q_P  = zeros(1, NumToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, NumToProcess);
trackResults.dllDiscrFilt   = inf(1, NumToProcess);
trackResults.pllDiscr       = inf(1, NumToProcess);
trackResults.pllDiscrFilt   = inf(1, NumToProcess);

% Remaining code and carrier phase for each tracking update
trackResults.remCodePhase   = inf(1, NumToProcess);
trackResults.remCarrPhase   = inf(1, NumToProcess);

% C/No and PLL lock detector of data channel
trackResults.DataCNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.DataPLD  = zeros(1,floor(NumToProcess/settings.CNoInterval));

% C/No and PLL lock detector of pilot channel
trackResults.PilotCNo = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.PilotPLD  = zeros(1,floor(NumToProcess/settings.CNoInterval));
trackResults.B1C_CNo  = zeros(1,floor(NumToProcess/settings.CNoInterval));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Consider the subcarrier of the BOC(1,1) modulation
codeLength = settings.codeLength;

% Summation interval
PDIcode = settings.intTime;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, ...
    1.0);

%--- PLL variables --------------------------------------------------------
% Calculate filter coefficient values
[pf3,pf2,pf1] = calcLoopCoefCarr(settings);

% -------- Number of acqusired signals ------------------------------------
TrackedNr =0 ;
for channelNr = 1:settings.numberOfChannels
    if channel(channelNr).status == 'T'
        TrackedNr = TrackedNr+1;
    end
end

% Start waitbar
hwb = waitbar(0,'Tracking...');

%Adjust the size of the waitbar to insert text
CNoPos=get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels

    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;

        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample)
        fseek(fid, ...
            dataAdaptCoeff*(settings.skipNumberOfBytes + channel(channelNr).codePhase-1), ...
            'bof');

        % Get a vector with the B1C data-channel code sampled 1x/chip
        B1CData = generateDataBOC11(settings,channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        B1CData = [B1CData(codeLength*2) B1CData B1CData(1)]; %#ok<AGROW>

        % Get a vector with the pilot BOC(1,1) spreading waveform
        pilotBOC11 = generatePilotBOC11(settings,channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        pilotBOC11 = [pilotBOC11(codeLength*2) pilotBOC11 pilotBOC11(1)];  %#ok<AGROW>

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO,
        codeFreq      = channel(channelNr).codeFreq;

        % Define residual code phase (in chips)
        remCodePhase  = 0.0;

        % Define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;

        % Define residual carrier phase
        remCarrPhase  = 0.0;

        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;

        % Carrier/Costas loop parameters
        d2CarrError  = 0.0;
        dCarrError   = 0.0;

        % For C/No computation
        CNoValue = zeros(1,3);
        tempCNoValue = zeros(1,3);

        %=== Process the number of specified code periods =================
        for loopCnt =  1:NumToProcess

            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 200 ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 20) == 0)

                Ln = newline;
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(TrackedNr),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt*10), ...
                    ' of ', int2str(NumToProcess*10), ' msec',Ln...
                    'Data C/No: ',int2str(CNoValue(1)),' (dB-Hz);',...
                    '   Pilot C/No: ',int2str(CNoValue(2)),' (dB-Hz)'];

                try
                    waitbar(loopCnt/NumToProcess, ...
                        hwb, ...
                        trackingStatus);
                catch %#ok<CTCH>
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

            %% Read next block of data ------------------------------------------------
            % Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) =(ftell(fid))/dataAdaptCoeff;

            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq;

            % Find the size of a "block" or code period in whole samples
            blksize = ceil((codeLength-remCodePhase) / codePhaseStep);

            % Read in the appropriate number of samples to process this
            % interation
            [rawSignal, samplesRead] = fread(fid, ...
                dataAdaptCoeff*blksize, settings.dataType);

            rawSignal = rawSignal';

            if (dataAdaptCoeff==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;  % transpose vector
            end

            % If did not read in enough samples, then could be out of
            % data - better exit
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                delete(hwb);
                return
            end

            %% Set up all the code phase tracking information -------------------------
            % Save remCodePhase for current correlation
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
            % Define index into early code vector
            tcode       = (remCodePhase-earlyLateSpc)*2 : ...
                codePhaseStep*2 : ...
                ((blksize-1) * codePhaseStep + remCodePhase- earlyLateSpc)*2;
            tcode2      = ceil(tcode) + 1;
            earlyCode   = B1CData(tcode2);

            % For pilot channel signal tracking
            p11_earlyCode   = pilotBOC11(tcode2);

            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc)*2 : ...
                codePhaseStep*2 : ...
                ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc)*2;
            tcode2      = ceil(tcode) + 1;
            lateCode    = B1CData(tcode2);

            % For pilot channel signal tracking
            p11_lateCode   = pilotBOC11(tcode2);

            % Define index into prompt code vector
            tcode       = remCodePhase*2 : ...
                codePhaseStep*2 : ...
                ((blksize-1)*codePhaseStep+remCodePhase)*2;
            tcode2      = ceil(tcode) + 1;
            promptCode  = B1CData(tcode2);

            % For pilot channel signal tracking
            p11_promptCode   = pilotBOC11(tcode2);

            remCodePhase = tcode(blksize)/2 + codePhaseStep - codeLength;

            %% Generate the carrier frequency to mix the signal to baseband -----------

            % Save remCarrPhase for current correlation
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;

            % Get the argument to sin/cos functions
            time    = (0:blksize) ./ settings.samplingFreq;
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));

            % Finally compute the signal to mix the collected data to
            % bandband
            carrsig = exp(-1i .* trigarg(1:blksize));

            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            iBasebandSignal = real(carrsig .* rawSignal);
            qBasebandSignal = imag(carrsig .* rawSignal);

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);

            % Correlation values for pilot BOC(1,1) spreading waveform
            p11_I_E = sum(p11_earlyCode  .* iBasebandSignal);
            p11_Q_E = sum(p11_earlyCode  .* qBasebandSignal);
            p11_I_P = sum(p11_promptCode .* iBasebandSignal);
            p11_Q_P = sum(p11_promptCode .* qBasebandSignal);
            p11_I_L = sum(p11_lateCode   .* iBasebandSignal);
            p11_Q_L = sum(p11_lateCode   .* qBasebandSignal);

            %% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);

            % Combined code tracking error estimation using data and pilot
            % chaannel signals
            % Pilot channel carrier phase is pi/2 rad ahead of the
            % data channel carrier phase.atan is not affectede by
            % the secondary code modulation
            p11_carrError = atan(-p11_I_P/p11_Q_P)/ (2.0 * pi);

            % Composite carrier tracking error
            carrError = (carrError *11 + p11_carrError*29)/40;

            % Implement carrier loop filter and generate NCO command
            d2CarrError = d2CarrError + carrError * pf3;
            dCarrError  = d2CarrError + carrError * pf2 + dCarrError;
            carrNco     = dCarrError + carrError * pf1;

            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;

            %% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_E ^2 + Q_E ^2) - sqrt(I_L ^2 + Q_L ^2)) / ...
                (sqrt(I_E ^2 + Q_E ^2) + sqrt(I_L ^2 + Q_L ^2)) * (1-earlyLateSpc);

            % Combined code tracking error estimation using data and pilot
            % chaannel signals
            p11_codeError = (sqrt(p11_I_E ^2 + p11_Q_E ^2) - sqrt(p11_I_L ^2 + p11_Q_L ^2)) / ...
                (sqrt(p11_I_E ^2 + p11_Q_E ^2) + sqrt(p11_I_L ^2 + p11_Q_L ^2))* (1-earlyLateSpc);
            % Composite code tracking error
            codeError = (codeError*11 + p11_codeError*29)/40;

            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;

            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;
            % Modify code freq based on NCO command
            codeFreq = channel(channelNr).codeFreq - codeNco;

            %% Record various measures to show in postprocessing ----------------------


            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;

            trackResults(channelNr).Pilot_I_P(loopCnt) = p11_I_P ;
            trackResults(channelNr).Pilot_Q_P(loopCnt) = p11_Q_P;
            %% CNo calculation --------------------------------------------------------

            if (rem(loopCnt,settings.CNoInterval)==0)
                % Computation of CNo and PLL detector output
                [CNoValue, PllDetector]= ...
                    Calc_CNo_PLD(trackResults(channelNr),settings,loopCnt);

                CNoCnt = loopCnt/settings.CNoInterval;

                % Save C/No for data channel: a o.5-0.5 filter is used to
                % smooth the results
                trackResults(channelNr).DataCNo(CNoCnt) = ...
                    CNoValue(1) * 0.5 + tempCNoValue(1) * 0.5;
                % Save PLL lock detector output for data channel
                trackResults(channelNr).DataPLD(CNoCnt) = PllDetector(1);

                % Save C/No and PLL lock detector output for pilot channel
                trackResults(channelNr).PilotCNo(CNoCnt) = ...
                    CNoValue(2) * 0.5 + tempCNoValue(2) * 0.5;
                trackResults(channelNr).B1C_CNo(CNoCnt) = ...
                    CNoValue(3) * 0.5 + tempCNoValue(3) * 0.5;
                trackResults(channelNr).PilotPLD(CNoCnt) = PllDetector(2);
            end
            tempCNoValue = CNoValue;

        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;

    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
