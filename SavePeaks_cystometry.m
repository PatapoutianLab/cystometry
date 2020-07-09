
% Script to find and plot bladder pressure peaks and sphincter EMG data
% Will save individual peak files (pressure and EMG) for later analysis
% Should be used before PlotandAnalyze file
%by Kara Marshall and Jason Keller
clear all; close all;

%Path to files
Path = '/Users/path/';

fontSz = 15;
sampleRate = 10000; %sampling rate per sec of matlab data
plotRaw = 1; % if 1, will plot bladder pressure and EMG
preSecToPlot = 20; %in seconds
postSecToPlot = 20; %in seconds
preSamplesToPlot = preSecToPlot*sampleRate; 
postSamplesToPlot = postSecToPlot*sampleRate; 
includeBladder = 1;
prePeakToAnalyze = 5; %in seconds
prePeakToAnalyze = 5; %in seconds
stimLength = 5;  %in seconds
downSampFactor = 5; %downsample to make plotting easier

cmapBladder = [0 0 0];
cmapSphincter = [0 0 0];

% mouseNums = {put each matlab file name, without extension, between apostrophes, with a space between them}
mouseNums = {'file_1'};
totalMice = size(mouseNums, 2);


emgData = {};
emgDataRms = {};
bladderData = {};


for i = 1:totalMice
    mouseNum = mouseNums{1,i};
    matName = [Path, mouseNum, '.mat'];
    load(matName);
   
    savename1 = [Path, mouseNum, '_peaksandlocs.xlsx'];
    savename2 = [Path, mouseNum, '_peakwidths.xlsx'];
    
    
    bladderData = data(:,1);
    emgData = data(:,2);
    emgDataRms = data(:,3);

    clear data;  %free up some memory
    minSampBetweenPeaks = sampleRate*20; %in seconds
    
    % find bladder pressure peaks above threshold:
    % first optionally clip data if artifacts at beginning/end
    timescaleRaw = 0:1/sampleRate:(size(bladderData,1)-1)/sampleRate;   % timescale in seconds
    startSec = 1; 
    startIndRaw = find(timescaleRaw > startSec, 1, 'first');
    endIndRaw = length(timescaleRaw);
    bladderData = bladderData(startIndRaw:endIndRaw) - min(bladderData(startIndRaw:endIndRaw)); %set minimum to zero, clip at beginning/end
    [peakLevels, bladderPeaks, w] = findpeaks(bladderData, 'MINPEAKHEIGHT', 15, 'MINPEAKDISTANCE', minSampBetweenPeaks); %change 'MINPEAKHEIGHT' to preferred threshold for a peak
 
    indToDelete1 = find(bladderPeaks < preSamplesToPlot); %get rid of timing too close to start for levels, to match the # of pks 
    indToDelete2 = find((length(bladderData) - bladderPeaks) < postSamplesToPlot);
    indToDelete = unique([indToDelete1; indToDelete2]);
    bladderPeaks(indToDelete)= []; %delete element if not long enough time period passed, or not enough time after
    peakLevels(indToDelete)= [];
    w(indToDelete) = [];
    
% Get locations (in seconds) and values of peaks to calculate the interpeak intervals and maximum pressures 
    peakLevelstable = array2table(peakLevels);%,'VariableNames', {'Peak location (s)'});
    peakWidthstable = array2table(w/sampleRate); % peak widths
    PeaksandLevels = [peakLevelstable bladderPeakstable];
    Peakwidths = [peakWidthstable]; % did not use this analysis
    
    writetable(PeaksandLevels, savename1);
    writetable(Peakwidths, savename2);
    
    
    indToDelete1 = find(bladderPeaks < preSamplesToPlot); %get rid of timing too close to start
    indToDelete2 = find((length(bladderData) - bladderPeaks) < postSamplesToPlot);
    indToDelete = unique([indToDelete1; indToDelete2]);
    bladderPeaks(indToDelete)= []; %delete element if not long enough time period passed, or not enough time after
    
 
    for k = 1:length(bladderPeaks) % make the single-peak files for later analysis and plotting
        startInd(k) = bladderPeaks(k) - preSamplesToPlot;
        endInd(k) = bladderPeaks(k) + postSamplesToPlot;
        emgTrace = emgData(startInd(k):endInd(k));
        bladderTrace = bladderData(startInd(k):endInd(k));
        timescaleSec = -preSecToPlot:1/sampleRate:postSecToPlot;
        savePulseName = [Path, mouseNum, '_', num2str(k), '.mat'];
        save(savePulseName, 'emgTrace', 'bladderTrace', 'timescaleSec');
    end
    
    if plotRaw
        timescale = downsample(0:1/sampleRate:(size(bladderData,1)-1)/sampleRate, downSampFactor);
%         startSec = 1; 
%         endSec = 850; %Can define selected start and end seconds if desired
        startInd = find(timescale > startSec, 1, 'first');
%         endInd = find(timescale > endSec, 1, 'first'); %if desired end is not full length
        endInd = length(timescale);
        figure;
        subplot(2,1,1);
        bladderDown = downsample(bladderData, downSampFactor);
        bladderToPlot = bladderDown(startInd:endInd)-min(bladderDown(startInd:endInd));
%         plot(timescale(startInd:endInd)-timescale(startInd), bladderToPlot, 'Color', cmapBladder, 'LineWidth', 0.25)
        plot(timescale(startInd:endInd), bladderToPlot, 'Color', cmapBladder, 'LineWidth', 0.25)
        xlabel('time (sec)', 'FontSize', fontSz);
        ylabel('bladder pressure (mmHg)', 'FontSize', fontSz);
        xLims = get(gca, 'XLim');
        axis([xLims(1) xLims(2) -10 50]); %change axis depending on data

        subplot(2,1,2);
        emgDown = downsample(emgData, downSampFactor);
        emgToPlot = emgDown(startInd:endInd)-mean(emgDown);
        plot(timescale(startInd:endInd), emgToPlot, 'Color', cmapSphincter, 'LineWidth', 0.25)
        xlabel('time (sec)', 'FontSize', fontSz);
        ylabel('sphincter EMG (mV)', 'FontSize', fontSz);
        xLims = get(gca, 'XLim');
        axis([xLims(1) xLims(2) -.5 .5]); % change axis depending on data.  Male and female axes will be very different for urethral sphincter EMG
        
        
        saveas(gcf, [Path, mouseNum, '.jpg'], 'jpg');
    
    end
    
end
%plot2svg([emgPath, 'EmgBladderPlot.svg'], gcf)