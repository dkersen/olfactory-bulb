% Calculates the average increase in spiking as a function of each MC's
% number of connected GCs

load('MCFeedbackResults.mat')
load('glomeruli.mat')
load('fullNetwork.mat')

R = sum(network,2);

numOdors = 5;
X_odor = [];
Y_odor = [];
X_nodor = [];
Y_nodor = [];
odorResp = [];
mitralNum = length(glomArray);

connectData = [];
connectDataOdor = [];
connectDataNodor = [];
        
        

for i = 1:numOdors
    odorMits = [];
    for j = 1:mitralNum
        if sum(ismember(glomArray(j),odorGloms{i})) > 0
            odorMits = [odorMits; j];
        end
    end
    rates_pre = sum(spikeTrainRecord{1,i}(:,1667:end),2);
    rates_post = sum(spikeTrainRecord{2,i}(:,1667:end),2);
    diff = rates_post - rates_pre;
    fb_receivers = find(feedbackRecord{i});
    odor_FB = intersect(fb_receivers, odorMits);
    nodor_FB = setdiff(fb_receivers, odorMits);
    
    X_odor = [X_odor; rates_pre(odor_FB)];
    Y_odor = [Y_odor; diff(odor_FB)];
    X_nodor = [X_nodor; rates_pre(nodor_FB)];
    Y_nodor = [Y_nodor; diff(nodor_FB)];
    
    connectData = [connectData; R(fb_receivers)];
    connectDataOdor = [connectDataOdor; R(odor_FB)];
    connectDataNodor = [connectDataNodor; R(nodor_FB)];

    odorResp = [odorResp; rates_pre(odorMits)];
end

X_g = [X_odor; X_nodor];
Y_g = [Y_odor; Y_nodor];


maximum = max(X_g);
rangeVector = 0:maximum;
meanChange = zeros(1,length(rangeVector));
SEMChange = zeros(1,length(rangeVector));
max_value = zeros(1,length(rangeVector));
min_value = zeros(1,length(rangeVector));

for i = 1:length(rangeVector)
    meanChange(i) = mean(Y_g(X_g==rangeVector(i)));
    SEMChange(i) = std(Y_g(X_g==rangeVector(i))) / sqrt(length(Y_g(X_g==rangeVector(i))));
    if ~isnan(meanChange(i))
        max_value(i) = max(Y_g(X_g==rangeVector(i)));
        min_value(i) = min(Y_g(X_g==rangeVector(i)));
    end
end



maximum = max(X_odor);
rangeVectorO = 0:maximum;
meanChangeO = zeros(1,length(rangeVectorO));
SEMChangeO = zeros(1,length(rangeVectorO));

for i = 1:length(rangeVectorO)
    meanChangeO(i) = mean(Y_odor(X_odor==rangeVectorO(i)));
    SEMChangeO(i) = std(Y_odor(X_odor==rangeVectorO(i))) / sqrt(length(Y_odor(X_odor==rangeVectorO(i))));
end



maximum = max(X_nodor);
rangeVectorN = 0:maximum;
meanChangeN = zeros(1,length(rangeVectorN));
SEMChangeN = zeros(1,length(rangeVectorN));

for i = 1:length(rangeVectorN)
    meanChangeN(i) = mean(Y_nodor(X_nodor==rangeVectorN(i)));
    SEMChangeN(i) = std(Y_nodor(X_nodor==rangeVectorN(i))) / sqrt(length(Y_nodor(X_nodor==rangeVectorN(i))));
end

maxY = max(Y_g);
yrange = 0:maxY;
meanR = zeros(1,length(yrange));
SEMR = zeros(1,length(yrange));

for i = 1:length(yrange)
    meanR(i) = mean(connectData(Y_g == yrange(i)));
    SEMR(i) = std(connectData(Y_g == yrange(i)))/sqrt(length(connectData(Y_g==yrange(i))));
end

maxYO = max(Y_odor);
yrangeO = 0:maxYO;
meanRO = zeros(1,length(yrangeO));
SEMRO = zeros(1,length(yrangeO));

for i = 1:length(yrangeO)
    meanRO(i) = mean(connectDataOdor(Y_odor == yrangeO(i)));
    SEMRO(i) = std(connectDataOdor(Y_odor == yrangeO(i)))/sqrt(length(connectDataOdor(Y_odor==yrangeO(i))));
end

maxYNO = max(Y_nodor);
yrangeNO = 0:maxYNO;
meanRNO = zeros(1,length(yrangeNO));
SEMRNO = zeros(1,length(yrangeNO));

for i = 1:length(yrangeO)
    meanRNO(i) = mean(connectDataNodor(Y_nodor == yrangeNO(i)));
    SEMRNO(i) = std(connectDataNodor(Y_nodor == yrangeNO(i)))/sqrt(length(connectDataNodor(Y_nodor==yrangeNO(i))));
end

errorbar(yrangeO, meanRO, SEMRO, 'LineStyle','none', 'marker','o')
hold on
errorbar(yrangeNO, meanRNO, SEMRNO, 'LineStyle','none', 'marker','o')