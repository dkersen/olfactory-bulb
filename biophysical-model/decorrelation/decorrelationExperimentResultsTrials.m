% Calculates correlations between odor firing rates as a function of the
% window size, for the with and without GC conditions

load('DCresultsTrial.mat')

%  windowLengths = [1000, 900, 800, 700, 600, 500, 400, 300, 200, 150, 100, 50];
windowLengths = [1000];
mitralNum = 3550;

[~,numOdors] = size(spikeTrainRecord);

time_length = length(spikeTrainRecord{1,1});
numComparisons = numOdors*(numOdors-1)/2;

traces = cell(numComparisons, length(windowLengths));
meanValue = zeros(numComparisons, length(windowLengths));

compCount = 0;
for t = 1:length(windowLengths)
    timespan = 1:round(windowLengths(t)/2):time_length-windowLengths(t);
    localTrace = zeros(1,length(timespan));
    compCount = compCount + 1;
    for m = 1:numOdors
        for n = m+1:numOdors
            localCount = 0;
            for i = 1:length(timespan)
                outputFreq1 = sum(spikeTrainRecord{1,m}(:, timespan(i):timespan(i)+windowLengths(t)),2);
                outputFreq2= sum(spikeTrainRecord{1,n}(:, timespan(i):timespan(i)+windowLengths(t)),2);
                localTrace(i) = corr(outputFreq1,outputFreq2);
                 if timespan(i) > 1666 && timespan(i)+windowLengths(t) < 10000
                    meanValue(compCount,t) = meanValue(compCount,t) + localTrace(i);
                    localCount = localCount + 1;
                end
            end
            traces{compCount,t} = localTrace;
            meanValue(compCount,t) = meanValue(compCount,t)/localCount;
        end
    end
end

traces_NG = cell(numComparisons, length(windowLengths));
meanValue_NG = zeros(numComparisons, length(windowLengths));

compCount = 0;
for t = 1:length(windowLengths)
    timespan = 1:round(windowLengths(t)/2):time_length-windowLengths(t);
    localTrace = zeros(1,length(timespan));
    compCount = compCount + 1;
    for m = 1:numOdors
        for n = m+1:numOdors
            localCount = 0;
            for i = 1:length(timespan)
                outputFreq1 = sum(spikeTrainRecord{2,m}(:, timespan(i):timespan(i)+windowLengths(t)),2);
                outputFreq2= sum(spikeTrainRecord{2,n}(:, timespan(i):timespan(i)+windowLengths(t)),2);
                localTrace(i) = corr(outputFreq1,outputFreq2);
                 if timespan(i) > 1666 && timespan(i)+windowLengths(t) < 10000
                    meanValue_NG(compCount,t) = meanValue_NG(compCount,t) + localTrace(i);
                    localCount = localCount + 1;
                end
            end
            traces_NG{compCount,t} = localTrace;
            meanValue_NG(compCount,t) = meanValue_NG(compCount,t)/localCount;
        end
    end
end

