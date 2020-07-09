
% Script to plot and analyze individual bladder pressure peak files and sphincter EMG activity
% Use after SavePeaks file
%By Kara Marshall and Jason Keller

clear all; close all;

Path = '/Users/path/';

emgThresh = 3; % in microV for removing the many very low values in EMG trace for proper mean and sum calculations.
fontSz = 15;
sampleRate = 10000;
wndsz = 300; %for smoothing, ex. 100 samples is 10ms @ 10kHz.  For calculating the window size for use during Gaussian
plotSeparate = 1; % will plot each peak
plotIntensityMap = 1; % will plot heatmaps
prepeakSecToPlot = 5; % in seconds, distance from file start to isolate peak values in the middle
postpeakSecToPlot = 5; % in seconds, these post/prepeak are for getting means and sums from only around the bladder contraction, which is roughly middle 10s. Chosen based on bladder plots.
preSecToPlot = 20; %in seconds
postSecToPlot = 20; %in seconds
baselinenoise = 2; %in seconds
baselinenoiseToPlot = baselinenoise*sampleRate;
prepeakSamplesToPlot = prepeakSecToPlot*sampleRate;
postpeakSamplesToPlot = postpeakSecToPlot*sampleRate;
preSamplesToPlot = preSecToPlot*sampleRate; 
postSamplesToPlot = postSecToPlot*sampleRate; 
totalSamples = preSamplesToPlot + postSamplesToPlot + 1;
% peakonlySamples = preSamplesToPlot - prepeakSamplesToPlot + postpeakSamplesToPlot + 1;
downSampFactor = 1; %downsample to make plotting easier - for some reason SVG plots not working if too large
peak = round(totalSamples/2);
peakonlySamplesstart = peak - prepeakSamplesToPlot;
peakonlySamplesend = peak  + postpeakSamplesToPlot;


cmapBladder = [0 0 0];
cmapSphincter = [0 0 0];

%% list data
%list each peak file, without extension, between apostrophes
PeakNames = {'peak1' 'peak2' 'peak3'};

%% decide data
% examplePooledData = [wt13275]; %describe data
Peaks1 = [PeakNames];
savename = [Path, 'Practice.xlsx'];
figureName = 'practice';
cmap1 = cmapBladder;
totalPeaks = size(Peaks1, 2);

%% compute
% for filtering/smooothing EMG:
gaussFilter = gausswin(wndsz); % value found by guess & check to define the appropriate filter size
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that it all adds up to 1 (apply filter accross matrix equivalently)
validSamples = totalSamples-wndsz+1;

emgData = zeros(totalSamples, totalPeaks);
emgDataAbs = zeros(totalSamples, totalPeaks);
emgSmoothtoMin = zeros(validSamples, totalPeaks);
emgSmoothtoMean = zeros(validSamples, totalPeaks);
emgSmoothNorm = zeros(validSamples, totalPeaks);
bladderData = zeros(totalSamples, totalPeaks);
bladderNorm = zeros(totalSamples, totalPeaks);

emgPauses = false(validSamples, totalPeaks); %samples with RMS value < mean noise
emgPauseLength = zeros(totalPeaks,1 );
emgBurstLength = zeros(totalPeaks,1 );
emgRmsMeanNoise = zeros(totalPeaks,1 );
emgRmsMeanStd = zeros(totalPeaks,1 );

%read data files:
for i = 1:totalPeaks
    pulseNum = Peaks1{1,i};
    matName = [Path, pulseNum, '.mat'];
    load(matName);
    
    %take bladder median for first two seconds of trace
    bladderTraceZero = median(bladderTrace(1:baselinenoiseToPlot));
    %subtract bladder mean before stim
    bladderData(:,i) = bladderTrace - bladderTraceZero;
    bladderNorm(:,i) = bladderData(:,i)./max(bladderData(:,i)); % normalize to peak
    bladderDataforArea(:,i) = bladderTrace - min(bladderTrace); %get trace zeroed out for Trapz calculation
    
    emgTraceZero = mean(emgTrace(1:preSamplesToPlot));
    emgData(:,i) = emgTrace - emgTraceZero; %normalize to mean before stim
    emgDataAbs(:,i) = abs(emgTrace);
    
    temp = emgTrace.^2; %take squared amplitude before smoothing
    emgSmoothTemp = sqrt(conv(temp, gaussFilter, 'valid')); %smoothing filter & sqrt
%     emgRmsMeanNoise(i) = mean(emgSmoothTemp(1:preSamplesToPlot)); %calculate mean before stim
%     emgRmsMeanStd(i) = std(emgSmoothTemp(1:preSamplesToPlot)); %calculate std before stim
%     emgSmooth(:,i) = emgSmoothTemp - emgRmsMeanNoise(i); %subract mean before stim
    emgRmsMeanNoise(i) = mean(emgSmoothTemp(1:baselinenoiseToPlot)); %calculate mean noise at beginning of trace
    emgSmoothtoMean(:,i) = emgSmoothTemp - emgRmsMeanNoise(i); %subract mean noise
    emgRmsMinNoise(i) = min(emgSmoothTemp(1:baselinenoiseToPlot)); %calculate min value at beginning of trace
    emgSmoothtoMin(:,i) = emgSmoothTemp - emgRmsMinNoise(i); %subract minimum noise value to get all positive values for RMS EMG analysis
    
end

xDown = downsample(timescaleSec(1:validSamples),downSampFactor);
[meanPause, stdPause, semPause] = grpstats(emgPauseLength,[],{'mean','std','sem'});

%% plot

if plotSeparate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:totalPeaks
        pulseNum = Peaks1{1,i};
        figure;
        subplot(2,1,1);
        plot(timescaleSec, bladderData(:,i), 'Color', cmapBladder, 'LineWidth', 1)
        xlabel('time (sec)', 'FontSize', fontSz);
        ylabel('bladder pressure (mmHg)', 'FontSize', fontSz);
%         axis tight;
        axis([-preSecToPlot postSecToPlot -10 30]); 
        subplot(2,1,2);
        hold on;

        plot(timescaleSec, emgData(:,i), 'Color', cmapSphincter, 'LineWidth', 1)
        xlabel('time (sec)', 'FontSize', fontSz);
        ylabel('sphincter EMG (mV)', 'FontSize', fontSz);
        xLims = get(gca, 'XLim');
        axis([-preSecToPlot postSecToPlot -.5 .5]); %for raw EMG
        set(gca, 'FontSize', fontSz);
        saveas(gcf, [Path, pulseNum, '.jpg'], 'jpg');
        
    end
end

if plotIntensityMap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hBladderHeatmap = figure;
    bladderDownsampled = downsample(bladderData(1:validSamples,:),downSampFactor);
%   bladderDownsampledpeaks = downsample(bladderData(1:peakonlySamples,:),downSampFactor);
    %bladderX = (1:1:400001)'; %create x-axis for bladder data for the trapz calculations
    bladderIntegrals = trapz(bladderDataforArea)/length(bladderDataforArea); % trapz integral per sample - sample is each bladder peak 40s interval (40K)
    bladderIntegralsforSort = trapz(bladderDownsampled)/length(bladderDownsampled); % trapz integral per sample - will have some values as negative, looks better than the absolute area values 
    bladderRMS = rms(bladderDownsampled); % positive area total
    [~,maxRowIndeces] = sort(bladderIntegralsforSort,'descend'); %sort by decreasing bladder pressure area under curve (AUC) using trapz
    bladderDataSorted = bladderDownsampled(:,maxRowIndeces); % use these data if you want a sorted heatmap, which looks cleaner
    heatmap(bladderDownsampled'); 
%     heatmap(bladderDataNorm(1:validSamples,:)')
    colormap bone % can change in the colormap feature after plotting
    caxis([-10 30]);
    timescale = -preSecToPlot: (postSecToPlot+preSecToPlot)/length(bladderDownsampled) :postSecToPlot;
    x1 = find(timescale>-20, 1, 'first');  %find additional indeces for x axis
    x2 = find(timescale>-10, 1, 'first');
    x3 = find(timescale>0, 1, 'first');
    x4 = find(timescale>10, 1, 'first');
    x5 = length(timescale)-1;
    set(gca,'XTickLabel', {'-20' '10' '0' '10' '20'});
    set(gca,'XTick', [x1 x2 x3 x4 x5]);
    xlabel('peak pressure time (s)', 'FontSize', fontSz);
    title('Bladder pressure (mmHg)', 'FontSize', fontSz);
    colorbar
    
    hEmgHeatmap = figure;
%     heatmap((emgDataSmooth(1:validSamples,:).*1000)')
    emgDownsampled = downsample(emgSmoothtoMean(1:validSamples,:).*1000,downSampFactor);
    xDown2 = downsample(timescaleSec(1:validSamples),downSampFactor);
%     [~,maxRowIndeces] = sort(mean(emgDownsampled)); %comment out to use bladder sort here as well, or sort by mean EMG
    emgDataSorted = emgDownsampled(:,maxRowIndeces);
    heatmap(emgDownsampled')
    colormap(1-gray);
    caxis([0 30]); 
    set(gca,'XTickLabel', {'-20' '10' '0' '10' '20'});
    set(gca,'XTick', [x1 x2 x3 x4 x5]);    
    xlabel('peak pressure time (s)', 'FontSize', fontSz);
    title('RMS EMG (uV)', 'FontSize', fontSz);
    colorbar
    
%   plot the bladder integrals and RMS EMG data in a table 
    meanRMSEMG = mean(emgSmoothtoMin(1:validSamples,:),1);
%     RMSEMGsums = sum(emgSmooth(1:validSamples,:),1)./40; % average RMS EMG voltage per sec
    RMSEMGsums = sum(emgSmoothtoMean(1:validSamples,:),1); % sum RMS EMG voltage 
    meanRMSEMGpeak = mean(emgSmoothtoMin(peakonlySamplesstart:peakonlySamplesend,:),1); % mean RMS EMG voltage during middle 10s of peak
    RMSEMGsumspeak = sum(emgSmoothtoMean(peakonlySamplesstart:peakonlySamplesend,:),1); % sum RMS EMG voltage during middle 10s of peak
    peakPressures = max(bladderDownsampled,[],1); % find max values of bladder pressure peak for each contraction
   
    emgSmoothtoMin2 = emgSmoothtoMin.*1000; % make new emg values for deleting those values below threshold, convert to microV
    emgToDelete = find(abs(emgSmoothtoMin2)<emgThresh); %find non-signal values in EMG
    emgSmoothtoMin2(emgToDelete) = NaN ; %get rid of all of the non-signal values for mean and sum calc
    
    meanCleanEMG = nanmean(emgSmoothtoMin2(peakonlySamplesstart:peakonlySamplesend,:),1); % calculate mean without NaN 
    sumCleanEMG = nansum(emgSmoothtoMin2(peakonlySamplesstart:peakonlySamplesend,:),1); % calculate sum without NaN 
    
%   write out quantification of bladder pressure and EMG (during peak times)and peak pressure amplitudes to excel file
    combinedBPRMS = [bladderIntegrals; bladderRMS; RMSEMGsumspeak; meanRMSEMGpeak; peakPressures; meanCleanEMG; sumCleanEMG];
    peakquantification = array2table(combinedBPRMS,'RowNames',{'bladderIntegrals', 'bladderRMS', 'RMSEMGsumspeak', 'meanRMSEMGpeak', 'peakPressures', 'meanCleanEMG', 'sumCleanEMG'});
    writetable(peakquantification, savename, 'WriteRowNames', 1);
    
end

