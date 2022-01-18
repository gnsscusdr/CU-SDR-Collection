function [trackResults, channel]= tracking(fid, channel, settings)
% Performs code and carrier tracking for all channels.
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
% (C) Developed for GPS L2C SDR by Yafeng Li, Nagaraj C. Shivaramaiah
% and Dennis M. Akos.
% Based on the original framework for GPS C/A SDR by Darius Plausinaitis,
% Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
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
% Freq of the CM code:
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
if settings.pilotTRKflag
    % Outputs from the correlators of pilot channel (In-phase):
    trackResults.Pilot_I_P            = zeros(1, NumToProcess);
    trackResults.Pilot_I_E            = zeros(1, NumToProcess);
    trackResults.Pilot_I_L            = zeros(1, NumToProcess);
    % Outputs from the correlators of pilot channel (Quadrature-phase):
    trackResults.Pilot_Q_E            = zeros(1, NumToProcess);
    trackResults.Pilot_Q_P            = zeros(1, NumToProcess);
    trackResults.Pilot_Q_L            = zeros(1, NumToProcess);
end

% Loop discriminators
trackResults.dllDiscr       = inf(1, NumToProcess);
trackResults.dllDiscrFilt   = inf(1, NumToProcess);
trackResults.pllDiscr       = inf(1, NumToProcess);
trackResults.pllDiscrFilt   = inf(1, NumToProcess);
% Remain code and carrier phase for position calculation
trackResults.remCodePhase   = inf(1, NumToProcess);
trackResults.remCarrPhase   = inf(1, NumToProcess);
% C/No
trackResults.CNo.VSMValue = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval/20));
trackResults.CNo.VSMIndex = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval/20));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================
%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips), mutiplying 2 is to onsider the
% return-zero CM and CL codes
earlyLateSpc = settings.dllCorrelatorSpacing *2 ;
% Consider the return-zero CL and CM codes
codeLength = settings.codeLength *2;
% Summation interval
PDIcode = settings.intTime;
% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, 1.0);

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
CNoPos = get(hwb,'Position');
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
        fseek(fid,dataAdaptCoeff*(settings.skipNumberOfBytes + ...
            channel(channelNr).codePhase-1), 'bof');
        
        % Get a vector with the CM code sampled 1x/chip
        cmCode = generateCMcode(channel(channelNr).PRN,settings);
        % Then make it possible to do early and late versions
        cmCode = [cmCode(codeLength) cmCode cmCode(1)];  %#ok<AGROW>
        
        if settings.pilotTRKflag
            % the CL code phase corresponding to channel(channelNr).codePhase
            CLCodePhase = channel(channelNr).CLCodePhase;
            % Get a vector with the CL code sampled 1x/chip
            CLCode = generateCLcode(channel(channelNr).PRN,settings);
            % Then make it possible to do early and late versions
            CLCode = [CLCode(settings.CLCodeLength*2) CLCode CLCode(1)]; %#ok<AGROW>
        end
        
        %--- Perform various initializations ------------------------------ 
        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis * 2;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;
        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        %carrier/Costas loop parameters
        d2CarrError  = 0.0;
        dCarrError   = 0.0;
        %C/No computation
        vsmCnt  = 0;CNo = 0;
        
        %=== Process the number of specified code periods =================
        for loopCnt =  1:NumToProcess
            
            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 200ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 10) == 0)
                
                Ln=newline;
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(TrackedNr),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt*20), ...
                    ' of ', int2str(NumToProcess*20), ' msec',Ln...
                    'C/No: ',CNo,' (dB-Hz)'];
                
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
            [rawSignal, samplesRead] = fread(fid, dataAdaptCoeff*blksize, ...
                settings.dataType);
            rawSignal = rawSignal';
            
            if (dataAdaptCoeff==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;  %transpose vector
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
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase/2;
            
            % Define index into early code vector
            tcode       = (remCodePhase - earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = cmCode(tcode2);
            
            % For pilot channel signal tracking
            if settings.pilotTRKflag
                earlyCodeCL = CLCode(tcode2+ codeLength *(CLCodePhase-1));
            end
            
            % Define index into late code vector
            tcode       = (remCodePhase + earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep + remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = cmCode(tcode2);
            
            % For pilot channel signal tracking
            if settings.pilotTRKflag
                lateCodeCL   = CLCode(tcode2+ codeLength *(CLCodePhase-1));
            end
            
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep + remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = cmCode(tcode2);
            
            % For pilot channel signal tracking
            if settings.pilotTRKflag
                promptCodeCL   = CLCode(tcode2+ codeLength *(CLCodePhase-1));
            end
            
            remCodePhase = (tcode(blksize) + codePhaseStep) - codeLength;
            
            %% Generate the carrier frequency to mix the signal to baseband -----------
            
            % Save remCarrPhase for current correlation
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;
            
            % Get the argument to sin/cos functions
            time    = (0:blksize) ./ settings.samplingFreq;
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            % Finally compute the signal to mix the collected data to
            % bandband
            carrsig = exp(1i .* trigarg(1:blksize));
            
            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);
            
            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            
            % For pilot channel signal tracking
            if settings.pilotTRKflag
                I_ECL = sum(earlyCodeCL  .* iBasebandSignal);
                Q_ECL = sum(earlyCodeCL  .* qBasebandSignal);
                I_PCL = sum(promptCodeCL .* iBasebandSignal);
                Q_PCL = sum(promptCodeCL .* qBasebandSignal);
                I_LCL = sum(lateCodeCL   .* iBasebandSignal);
                Q_LCL = sum(lateCodeCL   .* qBasebandSignal);
            end
            
            %% Find PLL error and update carrier NCO ----------------------
            
            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);
            
            % For pilot channel signal tracking
            if settings.pilotTRKflag
                % atan is not affectede by the NH code modulation
                carrErrorCL = atan(Q_PCL/I_PCL)/ (2.0 * pi);
                
                % As the L5I and L5Q power is the same, so a simple avergae
                % is calculated for the carrier phase error estimate
                carrError = (carrError + carrErrorCL)/2;
            end
            
            % Implement carrier loop filter and generate NCO command
            d2CarrError = d2CarrError + carrError * pf3;
            dCarrError  = d2CarrError + carrError * pf2 + dCarrError;
            carrNco     = dCarrError + carrError * pf1;
            
            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            
            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;
            
            %% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_E^2 + Q_E^2) - sqrt(I_L^2 + Q_L^2)) / ...
                (sqrt(I_E^2 + Q_E^2) + sqrt(I_L^2 + Q_L^2));
            
            % For pilot channel signal tracking
            if settings.pilotTRKflag
                codeErrorCL = (sqrt(I_ECL^2 + Q_ECL^2) - sqrt(I_LCL^2 + Q_LCL^2)) / ...
                    (sqrt(I_ECL^2 + Q_ECL^2) + sqrt(I_LCL^2 + Q_LCL^2));
                codeError = (codeError + codeErrorCL)/2;
                
                CLCodePhase = CLCodePhase + 1;
                if (CLCodePhase >= 76)
                    CLCodePhase = 1;
                end
            end
            
            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq/2;
            
            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis*2 - codeNco;
            
            %% Record various measures to show in postprocessing ----------------------
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError/2;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco/2;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;
            
            % Data channel
            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
            
            % Pilot Channel
            if settings.pilotTRKflag
                trackResults(channelNr).Pilot_I_E(loopCnt) = I_ECL;
                trackResults(channelNr).Pilot_I_P(loopCnt) = I_PCL;
                trackResults(channelNr).Pilot_I_L(loopCnt) = I_LCL;
                trackResults(channelNr).Pilot_Q_E(loopCnt) = Q_ECL;
                trackResults(channelNr).Pilot_Q_P(loopCnt) = Q_PCL;
                trackResults(channelNr).Pilot_Q_L(loopCnt) = Q_LCL;
            end
            
            
            if (rem(loopCnt,settings.CNo.VSMinterval)==0)
                vsmCnt = vsmCnt+1;
                CNoValue=CNoVSM(trackResults(channelNr).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                    trackResults(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
                trackResults(channelNr).CNo.VSMValue(vsmCnt)=CNoValue;
                trackResults(channelNr).CNo.VSMIndex(vsmCnt)=loopCnt;
                CNo=int2str(CNoValue);
            end
            
        end % for loopCnt
        
        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;
        
    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
