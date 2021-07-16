% Calculates average correlations for same and different network external
% feedback conditions

load('GCFeedbackResults.mat')
load('glomeruli.mat')

[numConditions ,numTrials] = size(spikeTrainRecord);
record = cell(numConditions, numTrials);
for i = 1:numConditions
    for j = 1:numTrials
        record{i,j} = spikeTrainRecord{i,j}(:,1668:end);
    end
end

windowLengths = [100];
numWindowLengths = length(windowLengths);

indices = [];
for i = 1:length(odorGloms)
    indices = [indices, find(odorGloms(i) == glomArray)];
end

[~,time_length] = size(record{1,1});

coor_same = zeros(numTrials,numWindowLengths);
coor_diff = zeros(numTrials,numWindowLengths);
strength_odor = zeros(numTrials, numWindowLengths);
strength_diff = zeros(numTrials,numWindowLengths);
results_odor = zeros(numTrials,numWindowLengths);
results = zeros(numTrials,2);

for t = 1:numWindowLengths
    timespan = 1:round(windowLengths(t)/2):time_length-windowLengths(t);

    for k = 1:numTrials
        localTraceSame = zeros(1,length(timespan));
        localTraceDiff = zeros(1,length(timespan));
        for i = 1:length(timespan)
            Z1 = sum(record{1,k}(indices,timespan(i):timespan(i)+windowLengths(t)),2);
            Z2 = sum(record{2,k}(indices,timespan(i):timespan(i)+windowLengths(t)),2);
            Z3 = sum(record{3,k}(indices,timespan(i):timespan(i)+windowLengths(t)),2);
            Z4 = sum(record{4,k}(indices,timespan(i):timespan(i)+windowLengths(t)),2);
            Z5 = sum(record{5,k}(indices,timespan(i):timespan(i)+windowLengths(t)),2);
            localTraceSame(i) = corr(Z2-Z1,Z3-Z1);
            localTraceDiff(i) = corr(Z2-Z1,Z5-Z4); 
        end
        coor_same(k,t) = mean(localTraceSame);
        coor_diff(k,t) = mean(localTraceDiff);
    end
end

