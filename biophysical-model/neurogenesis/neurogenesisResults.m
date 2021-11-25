% Calculates the correlation changes for each round of neurogenesis

% load('controlResults.mat')

[numGeneses, numOdors] = size(spikeTrainRecord);

%numGeneses = 6;
[mitralNum, ~] = size(odorAmpMatrix);
[~, time] = size(spikeTrainRecord{1,1});

windowLength = 50;
time_length = length(spikeTrainRecord{1,1});

numComparisons = numOdors*(numOdors-1)/2;
odorCompIdentity = zeros(numComparisons,2);
compCount = 0;
for i = 1:numOdors
    for j = i+1:numOdors
        compCount = compCount + 1;
        odorCompIdentity(compCount,:) = [i,j];
    end
end
        
meanValue = zeros(numGeneses, numComparisons);
meanChange = zeros(numGeneses-1,numComparisons);

for g = 1:numGeneses
    compCount = 0;
    timespan = 1:round(windowLength/2):time_length-windowLength;
    for m = 1:numOdors
        for n = m+1:numOdors
            localTrace = zeros(1,length(timespan));
            compCount = compCount + 1;
            localCount = 0;
            for i = 1:length(timespan)
                outputFreq1 = sum(spikeTrainRecord{g,m}(:, timespan(i):timespan(i)+windowLength),2);
                outputFreq2= sum(spikeTrainRecord{g,n}(:, timespan(i):timespan(i)+windowLength),2);
                x = corr(outputFreq1,outputFreq2);
                if ~isnan(x)
                    localTrace(i) = x;
                end
                 if timespan(i) > 1667 && timespan(i)+windowLength < 3333 %looking at the second sniff
                    meanValue(g,compCount) = meanValue(g,compCount) + localTrace(i);
                    localCount = localCount + 1;
                end
            end
            meanValue(g,compCount) = meanValue(g,compCount)/localCount;
            if g > 1
                meanChange(g-1, compCount) = meanValue(g,compCount) - meanValue(1,compCount);
            end
        end
    end
end



