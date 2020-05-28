function acqResults = acquisition(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for B2a signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - (=20+X+1) ms of raw signal from the 
%                       front-end.The first 20+X ms segment is in order
%                       to include at least the first Xms of a CM code; 
%                       The last 1ms is to make sure the index does not
%                       exceeds matrix dimensions of 10ms long.
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Developed for BDS B2a SDR by Yafeng Li, Nagaraj C. Shivaramaiah 
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
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $


%% Condition input signal to speed up acquisition ===========================

% If input IF signal freq. is too high, a resampling strategy is applied
% to speed up the acquisition. This is user selectable.
if (settings.samplingFreq > settings.resamplingThreshold && ...
                                       settings.resamplingflag == 1)
   
    %--- Filiter out signal power outside the main lobe of CM code ------------------
    fs = settings.samplingFreq;
    IF = settings.IF;
    % Bandwidth of CM mian lobe: 0.5e6 is a margin to make sure most of CA
    % code power will not filtered out
    BW = settings.codeFreqBasis*2 + 0.5e6;
    % Filter parameter
    w1 = (IF)-BW/2;
    w2 = (IF)+BW/2;
    wp = [w1*2/fs-0.002 w2*2/fs+0.002];
    % Filter coefficients
    b  = fir1(700,wp);
    % Filter operation
    longSignal = filtfilt(b,1,longSignal);
    
    % --- Find resample frequency -----------------------------------------
    % Refer to bandpass sampling theorem(Yi-Ran Sun,Generalized Bandpass
    % Sampling Receivers for Software Defined Radio)
    
    % Upper boundary frequency of the bandpass IF signal
    fu = settings.IF + BW/2;
    
    % Lower freq. of the acceptable sampling Freq. range
    n = floor(fu/BW);
    if (n<1)
        n = 1;
    end
    lowerFreq = 2*fu/n;
    
    % Lower boundary frequency of the bandpass IF signal
    fl = settings.IF - BW/2;
    
    % Upper boundary frequency of the acceptable sampling Freq. range
    if(n>1)
        upperFreq = 2*fl/(n-1);
    else
        upperFreq = lowerFreq;
    end
    
    % Save orignal Freq. for later use
    oldFreq = settings.samplingFreq;
    
    % Take the center of the acceptable sampling Freq. range as
    % resampling frequency. As settings are used to generate local
    % B2a code samples, so assign the resampling freq. to settings.
    % This can not change the settings.samplingFreq outside this 
    % acquisition function.
    settings.samplingFreq = ceil((lowerFreq + upperFreq)/2);
    
    %--- Downsample input IF signal -------------------------------------------------
    % Signal length after resampling
    signalLen = floor((length(longSignal)-1) /oldFreq * settings.samplingFreq);
    % Resampled signal index
    index = ceil((0:signalLen-1)/settings.samplingFreq * oldFreq);
    index(1) = 1;
    % Resampled signal
    longSignal = longSignal(index);
    
    % For later use
    oldIF = settings.IF;
    
    % Resampling is equivalent to down-converting the original IF by integer
    % times of resampling freq.. So the IF after resampling is equivalent to:
    settings.IF = rem(settings.IF,settings.samplingFreq);

end % resampling input IF signals

%% Acquisition initialization =============================================

%--- Find number of samples for differnet long signals ------------------------------
% Number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));              

% Number of samples of 2 ms L5 spreading codes
len2ms = samplesPerCode*2;

% Find number of samples per CM code chip
samples2CodeChip   = ceil(settings.samplingFreq / settings.codeFreqBasis)*2;

%--- Cut 2 cm input signal to do correlate ----------------------------------
sig2ms = longSignal(1:len2ms);

%--- Generate input and local signal to to correlate ------------------------
% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (len2ms-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band
numberOfFrqBins = round(settings.acqSearchBand * 2 / settings.acqStep) + 1;

%--- Initialize arrays to speed up the code -----------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, len2ms);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);

%--- Initialize acqResults and related variables -----------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, max(settings.acqSatelliteList));
% PRN code phases of detected signals
acqResults.codePhase    = zeros(1, max(settings.acqSatelliteList));
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, max(settings.acqSatelliteList));

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList

%% Correlate signals ======================================================   
    
    % Generate all L5I + L5Q code and sample them according to the sampling freq.
    B2aDataTable  = makeB2aDataTable(PRN,settings);
    B2aPilotTable = makeB2aPilotTable(PRN,settings);
    
    % generate local code duplicate to do correlate
    localB2aDataCode = [B2aDataTable, zeros(1,samplesPerCode)];
    localB2aPilotCode = [B2aPilotTable, zeros(1,samplesPerCode)];
    
    %--- Perform DFT of PRN code ------------------------------------------
    B2aDataFreqDom = conj(fft(localB2aDataCode));
    B2aPilotFreqDom = conj(fft(localB2aPilotCode));
    
    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins

        %--- Generate carrier wave frequency grid  -----------------------
        frqBins(frqBinIndex) = settings.IF - settings.acqSearchBand + ...
                               settings.acqStep * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(1i*frqBins(frqBinIndex) * phasePoints);
        
        %--- "Remove carrier" from the signal -----------------------------
        I1      = real(sigCarr .* sig2ms);
        Q1      = imag(sigCarr .* sig2ms);

        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom1 = fft(I1 + 1i*Q1);

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        convB2aData = IQfreqDom1 .* B2aDataFreqDom;
        convB2aPilot = IQfreqDom1 .* B2aPilotFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        results(frqBinIndex, :) = abs(ifft(convB2aData)) + abs(ifft(convB2aPilot));
    
    end % frqBinIndex = 1:numberOfFrqBins
  
%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [~, frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize, codePhase] = max(max(results));

    %--- Find 1 chip wide B2a code phase exclude range around the peak ----
    excludeRangeIndex1 = codePhase - samples2CodeChip;
    excludeRangeIndex2 = codePhase + samples2CodeChip;
    % Exclude range around the peak within X ms CM code period.
    excludeRangeIndex3 = codePhase - samplesPerCode + samples2CodeChip;
    excludeRangeIndex4 = codePhase + samplesPerCode - samples2CodeChip;

    %--- Correct PRN code phase exclude range if the range includes array
    %boundaries
    % The left including range before the peak
    if excludeRangeIndex1 < 1
        leftRange = [];        
    else
        leftRange = max(1,excludeRangeIndex3): excludeRangeIndex1;
    end
    % The right including range after the peak
    if excludeRangeIndex2 >= len2ms
        rightRange = [];
    else
        rightRange = excludeRangeIndex2:min(excludeRangeIndex4,len2ms);
    end
    
    % Combine the left and right Ranges together
    codePhaseRange = [leftRange rightRange];

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold
%% Fine resolution frequency search =======================================

        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Set related parameters and initialize arrays to speed up the code ------
        
        % Number of the frequency bins for fine acquisition: use 400Hz fine
        % acquisition band, and 25Hz step
        NumOfFineBins = round(settings.acqStep / 25) + 1;
        
        % Carrier frequencies of the frequency bins
        FineFrqBins     = zeros(1, NumOfFineBins);
        
        % Search results of all frequency bins
        FineResult = zeros(1,NumOfFineBins);
        
        %--- Generate B2a data and pilot codes sequence -------------------
        % settings.fineNoncoh ms long 
        B2aDataCode = generateB2aDataCode(PRN,settings);
        B2aPilotCode = generateB2aPilotCode(PRN,settings);
        
        % Sampling index
        codeValueIndex = floor((ts * (1:settings.fineNoncoh*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
        
        % Sampled data and pilot codes
        longB2aDataCode = B2aDataCode((rem(codeValueIndex, settings.codeLength) + 1));  
        longB2aPilotCode = B2aPilotCode((rem(codeValueIndex, settings.codeLength) + 1));
        
        %--- Find phase points of the local carrier wave -------------------
        finePhasePoints = (0 : (settings.fineNoncoh*samplesPerCode-1)) * 2 * pi * ts;
        
        % 10cm incoming signal
        sigFineACQ = longSignal(codePhase:codePhase + settings.fineNoncoh*samplesPerCode -1);
        
        % Coherent integration for each code
        sumPerCode1 = zeros(1,settings.fineNoncoh);
        sumPerCode2 = zeros(1,settings.fineNoncoh);
                
        %--- Search different frequency bins -------------------------------
        for FineBinIndex = 1 : NumOfFineBins
            
            % Carrier frequencies of the frequency bins
            FineFrqBins(FineBinIndex) = frqBins(frequencyBinIndex) -...
                             settings.acqStep/2 + 25 * (FineBinIndex - 1);
            % Generate local sine and cosine
            sigCarr20cm = exp(1i*FineFrqBins(FineBinIndex) * finePhasePoints);
            
            % Wipe off B2a code and carrier from incoming signals to 
            % produce baseband signal. This is for data channel
            basebandSig1 = longB2aDataCode .* sigCarr20cm .* sigFineACQ;
            
            % This is for pilot channel
            basebandSig2 = longB2aPilotCode .* sigCarr20cm .* sigFineACQ;
            
            % Non-coherent integration for each code
            for index = 1:settings.fineNoncoh
                sumPerCode1(index) = sum( basebandSig1( samplesPerCode * ...
                                     (index-1)+1:samplesPerCode*index ) );
                sumPerCode2(index) = sum( basebandSig2( samplesPerCode * ...
                                     (index-1)+1:samplesPerCode*index ) );                  
            end
            
            FineResult(FineBinIndex) = sum(abs(sumPerCode1)) + sum(abs(sumPerCode2));

        end % FineBinIndex = 1 : NumOfFineBins
        
        % Find the fine carrier freq. -------------------------------------
        % Corresponding to the largest noncoherent power
        [~,maxFinBin] = max(FineResult);
        acqResults.carrFreq(PRN) = FineFrqBins(maxFinBin);
        
        % Code phase acquisition result
        acqResults.codePhase(PRN) = codePhase;
        
        %signal found, if IF = 0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end

%% Downsampling recovery ============================================
       % Find acquisition results corresponding to orignal sampling freq
        if (exist('oldFreq', 'var') && settings.resamplingflag == 1)
            % Find code phase
            acqResults.codePhase(PRN) = floor((codePhase - 1)/ ...
                                        settings.samplingFreq *oldFreq)+1;
            
            % Doppler frequency
            if (settings.IF >= settings.samplingFreq/2)
                % In this condition, the FFT computed freq. is symmetric
                % with the true frequemcy about half of the sampling
                % frequency, so we have the following:
                IF_temp = settings.samplingFreq - settings.IF;
                doppler = IF_temp - acqResults.carrFreq(PRN);
            else
                doppler = acqResults.carrFreq(PRN) - settings.IF;
            end
            % Carrier freq. corresponding to orignal sampling freq
            acqResults.carrFreq(PRN) = doppler + oldIF;
        end
           
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
