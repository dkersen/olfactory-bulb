% Explores synaptic parameters to tune results to lateral inhibition
% experiments from Arevian et al (2008)
rng shuffle
% 
network = load('fullNetworkSlice.mat');
network = network.('network');
sliceDistance = load('distanceSlice.mat');
sliceDistance = sliceDistance.('distance'); 


mitPairs = [1493,2387;570,2023;931,1575;755,770;537,1803;570,959;1010,2250;879,1283;537,1763;863,1223;135,2250;2248,2250;833,2103;856,865;761,2250;931,2231];
[mitralNum, granuleNum] = size(network);


mParam = zeros(8, mitralNum);
gParam = zeros(8, granuleNum);

% Establish the intrinsic mitral cell parameters
for i = 1:mitralNum
    while mParam(1,i) >= mParam(2,i)
        mParam(1,i) = -58; % mVr (resting potential)
        mParam(2,i) = -49; % mVt (threshold potential)
    end
end
mParam(3,:) = 0.02*ones(1,mitralNum); %a (Izhikevich parameter a)
mParam(4,:) = 12*ones(1,mitralNum); %b_m (Izhikevich parameter b)
mParam(5,:) = 2.5* ones(1,mitralNum); %k_m (Izhikevich parameter k)

for i = 1:mitralNum
    while mParam(6,i) >= mParam(2,i)
        mParam(6,i) = -65; % c_m (Izhikevich parameter c)
    end
end
mParam(7,:) = 13*ones(1,mitralNum); %d_m (Izhikevich parameter d)  
mParam(8,:) = 191*ones(1,mitralNum); % cap_m (capacitance)
%mParam(9,:) = normrnd(0.00113,0.00113/10,1,mitralNum);




%recurrent MC synaptic parameters 
tauG_m = 18*ones(mitralNum, granuleNum); %tau_GABA

%odor input MC parameters (external)
mAMPA = normrnd(6.7,6.7/10,mitralNum, 1); % gAMPAx 
tauA_m = normrnd(14.3108,14.3108/10,mitralNum, 1); %tauAMPAx 
mNMDA = normrnd(12,12/10,mitralNum, 1); %gNMDAx
tauNr_m = normrnd(13,13/10,mitralNum, 1); %tauNrx
tauNd_m = normrnd(70,70/10,mitralNum, 1); %tauNdx

% Establish the intrinsic granule cell parameters
for i = 1:granuleNum
    while gParam(1,i) >= gParam(2,i) || gParam(2,i) - gParam(1,i) < 20
        gParam(1,i) = -71; %gVr (resting potential)
        gParam(2,i) = -39; %gVt (threshold potential)
    end
end
gParam(3,:) = 0.01*ones(1,granuleNum); %a_g (Izhikevich parameter a)

%Establish Izhikevich parameters based on rheobase and input resistance
for i = 1:granuleNum
    rheobase = 0;
    inres = -1;
    while (rheobase < 10 || rheobase > 70 || inres < 0.25 || inres > 1.5) || b > 0
        b = -2/15; 
        k = 1/15;
        rheobase = (b+k*(-gParam(1,i)+gParam(2,i)))^2/(4*k);
        inres = 1/(b - k*(gParam(1,i)-gParam(2,i)));
    end
    gParam(4,i) = b; %b_g (Izhekevich parameter b)
    gParam(5,i) = k; %k_g (Izhikevich parameter k)
end
for i = 1:granuleNum
    while gParam(6,i) >= gParam(2,i)
        gParam(6,i) = -75; %c_g (Izhikevich parameter c)
    end
end

gParam(7,:) = 1.2*ones(1,granuleNum); %d_g (Izhikevich parameter d)
gParam(8,:) = 48*ones(1,granuleNum); %cap_g (capacitance)

%Recurrent GC feedback parameters
tauA_g = 5.5*ones(mitralNum,granuleNum); %tau_AMPA
tauNr_g = 10*ones(mitralNum,granuleNum); %tau_NMDA_rise
tauNd_g = 80*ones(mitralNum,granuleNum); %tau_NMDA_decay

[len,~] = size(mitPairs);
networks = cell(1,len);
distances = cell(1,len);
mParameters = cell(3,len);
gParameters = cell(4,len);
removals = cell(1,len);

for i = 1:len
    remove = [];
    networks{i} = [network(mitPairs(i,1),:); network(mitPairs(i,2),:)];
    
    mParameters{1,i} = mParam(:, mitPairs(i,:));
    %Blank for x
    mParameters{3,i} = [tauG_m(mitPairs(1,1),:);tauG_m(mitPairs(i,2),:)];    

    gParameters{1,i} = gParam;
    %Blank for x
    gParameters{2,i} = [tauA_g(mitPairs(i,1),:); tauA_g(mitPairs(i,2),:)];
    gParameters{3,i} = [tauNr_g(mitPairs(i,1),:); tauNr_g(mitPairs(i,2),:)];
    gParameters{4,i} = [tauNd_g(mitPairs(i,1),:); tauNd_g(mitPairs(i,2),:)];

    count = 1;
    for j = 1:granuleNum
        if sum(networks{i}(:,j))==0
            remove(count) = j;
            count = count + 1;
        end
    end
    networks{i}(:,remove) = [];
    mParameters{3,i}(:,remove) = [];
    
    gParameters{1,i}(:,remove) = [];
    gParameters{2,i}(:,remove) = [];
    gParameters{3,i}(:,remove) = [];
    gParameters{4,i}(:,remove) = [];
    removals{i} = remove;
end

f = @(x)networkCost(x,network, sliceDistance, mitPairs, networks,mParameters, gParameters, removals);

numTrials = 10000;


for i = 1:numTrials
    
% 
    x0(1) = 0.73;
    x0(2) = 0.84;
    x0(3) = 0.13;
    x0(4) = 0.006;
    x0(5) = 675;



%     x0(1) = rand;
%     x0(2) = rand;
%     x0(3) = rand;
%     x0(4) = rand*0.1;
%     x0(5) = 675;

    cost = networkCost(x0,network, sliceDistance, mitPairs, networks, mParameters, gParameters, removals);
    disp(x0)
	disp(cost)
     if i == 1 || cost < error
        error = cost;
        xopt = x0;
        fname = sprintf('result%d.mat',i);
        save(fname,'xopt','cost');
    end
    disp('I am min');
    disp(xopt);
    disp(error);
end





function cost = networkCost(x,network, paraDistance, mitPairs, networks, mParameters, gParameters, removals)
    
    [mitralNum, granuleNum] = size(network);
    L = x(5)*ones(mitralNum,1);
    Lmat = repmat(L,[1,granuleNum]);
    mGABA = x(1) * exp(-paraDistance./Lmat);
    gAMPA = x(2);
    gNMDA = x(3);

    [len,~] = size(mitPairs);

    for i = 1:len
        mParameters{2,i} = [mGABA(mitPairs(1,1),:); mGABA(mitPairs(i,2),:)];    
        mParameters{2,i}(:,removals{i}) = [];
    end
    
    cost = sliceCost(x, mitPairs, networks, mParameters, gParameters, gAMPA, gNMDA); 
end

function sCost = sliceCost(x, mitPairs, networks, mParameters, gParameters, gAMPA, gNMDA)
    [len,~] = size(mitPairs);
    costArray = zeros(1,len);
%     inputArray = [700 700  680 680 680
%                   130 130  130 130 120 
%                   1280 1300 1250 1250 1250];
    for i = 1:len
        mParam = mParameters{1,i};
        mGABA = mParameters{2,i};
        tauG_m = mParameters{3,i};

        gParam = gParameters{1,i};
        tauA_g = gParameters{2,i};
        tauNr_g = gParameters{3,i};
        tauNd_g = gParameters{4,i};

        freq1 = latInhib(x,700,750,networks{i}, 1, mParam, mGABA, tauG_m, gParam, tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA);
        freq2 = latInhib(x,700,750,networks{i}, 2, mParam, mGABA, tauG_m, gParam, tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA);
        frac = 1-poslin(-0.00427839*freq1^2 + 0.590149 * freq1 - 6.06278)/100;
        costArray(i) = costArray(i) + (frac*freq1 - freq2)^2;
        if freq1 > 90
             costArray(i) = costArray(i) + (freq1 - 90)^2;
        elseif freq1 < 65
             costArray(i) = costArray(i) + (freq1 - 65)^2;
        end
        freq1 = latInhib(x,130,750,networks{i}, 1, mParam, mGABA, tauG_m, gParam, tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA);
        freq2 = latInhib(x,130,750,networks{i}, 2, mParam, mGABA, tauG_m, gParam, tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA);
        frac = 1-poslin(-0.00427839*freq1^2 + 0.590149 * freq1 - 6.06278)/100;
         costArray(i) = costArray(i) + (frac*freq1 - freq2)^2;
       if freq1 > 15
             costArray(i) = costArray(i) + (freq1 - 15)^2;
        elseif freq1 < 5
             costArray(i) = costArray(i) + (freq1 - 5)^2;
        end
        freq1 = latInhib(x,1500,750,networks{i}, 1, mParam, mGABA, tauG_m, gParam, tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA);
        freq2 = latInhib(x,1500,750,networks{i}, 2, mParam, mGABA, tauG_m, gParam, tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA);
        frac = 1-poslin(-0.00427839*freq1^2 + 0.590149 * freq1 - 6.06278)/100;
        costArray(i) = costArray(i) + (frac*freq1 - freq2)^2;   
       if freq1 > 130
             costArray(i) = costArray(i) + (freq1 - 130)^2;
        elseif freq1 < 100
             costArray(i) = costArray(i) + (freq1 - 100)^2;
        end
    end
    sCost = mean(costArray);
end

function f1 = latInhib(x, I1, I2, network, numFlag, mParam, mGABA, tauG_m, gParam, tauA_g, tauNr_g, tauNd_g, gAMPA, gNMDA)
    TS = 0.1;
    delayTime = 100;
    recordTime = 400;
    tspan = 0:TS:(delayTime+recordTime);

    [mitralNum,granuleNum] = size(network);
   

    k = x(4);
    W = 0.5;

    mVr = mParam(1,:); mVt = mParam(2,:);
    a_m = mParam(3,:); b_m = mParam(4,:); k_m = mParam(5,:); c_m = mParam(6,:); d_m = mParam(7,:);
    Cm = mParam(8,:);  
    
    gVr = gParam(1,:); gVt = gParam(2,:);
    a_g = gParam(3,:); b_g = gParam(4,:); k_g = gParam(5,:); c_g = gParam(6,:); d_g = gParam(7,:);
    Cg = gParam(8,:); 
    
    mV = ones(1,mitralNum) .* mVr;
    mU = zeros(1,mitralNum);
    gV = ones(1,granuleNum) .* gVr;
    gU = zeros(1,granuleNum);
    
    f1 = 0;
    mI1 = 0;
    mI2 = 0;
    
    if numFlag == 1
        mI1 = I1;
    elseif numFlag == 2
        mI1 = I1;
        mI2 = I2;
    end
        
    mIntegG = zeros(mitralNum, granuleNum);
    gIntegN = zeros(mitralNum, granuleNum);
    gIntegA = zeros(mitralNum, granuleNum);
    gX = zeros(mitralNum,granuleNum);
    
    
    numRecN = 1;
    numRecA = 1;


    for t = 1:length(tspan)
        mI(1) = mI1;
        mI(2) = mI2;
        mFired = find(mV >= 30);
        gFired = find(gV >= 25); 
        
        if mV(1) >= 30 && t > delayTime/TS
            f1 = f1 + 1/((recordTime)/1000);
        end
        
        gIntegA(mFired,:) = gIntegA(mFired,:) + network(mFired,:)*W.*(1-gIntegA(mFired,:));
        gX(mFired,:) = gX(mFired,:) + network(mFired,:)*W.*(1-gX(mFired,:));
        for i = 1:length(mFired)
            gRel = find(network(mFired(i),:));
            mIntegG(:,gRel) = mIntegG(:,gRel) + k*network(:,gRel)*W.*(1-mIntegG(:,gRel));
        end
        mIntegG(:,gFired) = mIntegG(:,gFired) + network(:,gFired)*W.*(1-mIntegG(:,gFired));
        
        Ei = -70; 
        Ee = 0;
        alphaM = 0.03;
        alphaG = 0.1;
        Mg = 0.2;


        %Remember parentheses
        Ig = mI-sum(mGABA .* mIntegG,2)' .* (mV - Ei);
        mI = Ig; % .* exp(-((Ia+In)-mu).^2/sigma);
        Ja = -(gV-Ee) .* (sum(gAMPA.*gIntegA));
        Jn = -(gV-Ee) .* ((sum(gNMDA.*gIntegN)).* 1./(1+exp(-0.062*gV)*Mg/3.57));
        gI = Ja + Jn;         
        
        mIntegG = mIntegG ./ exp(TS./tauG_m);
        
        gIntegA = gIntegA ./ exp(TS./tauA_g);
        gIntegN = gIntegN + TS*(-gIntegN./tauNd_g + alphaG*gX.*(1-gIntegN));
        gX = gX ./ exp(TS./tauNr_g);                
        
       [mV, mU, gV, gU] = izhikevich(mV, mU, gV, gU, mI, gI, mFired, gFired, TS,...
                a_m, b_m, c_m, d_m, k_m, a_g, b_g, c_g, d_g, k_g, mVr, mVt, gVr, gVt, Cm, Cg);   
    end
end

function [mI, gI, mIG,mIA,mIN,Xm,gIA,gIN, Xg] = current(mV, gV,...
    mIntegG, mIntegA, mIntegN, mX, gIntegA, gIntegN, gX,...
    Mg, mGABA, tauG_m, mAMPA, tauA_m, mNMDA, tauNr_m, tauNd_m,...
        gAMPA, tauA_g, gNMDA, tauNr_g, tauNd_g, TS)
        Ei = -70; 
        Ee = 0;
        alphaM = 0.03;
        alphaG = 0.1;
%         mu = 20;
%         sigma = 500;

        %Remember parentheses
        Ig = -sum(mGABA .* mIntegG,2)' .* (mV - Ei);
        Ia = - (mV-Ee) .* sum(mAMPA.*mIntegA,2)';
        In = -(mV-Ee) .* (sum(mNMDA.*mIntegN,2)') .* 1./(1+exp(-0.062*mV)*Mg/3.57);
        mI = Ig + Ia + In; % .* exp(-((Ia+In)-mu).^2/sigma);
        Ja = -(gV-Ee) .* (sum(gAMPA.*gIntegA));
        Jn = -(gV-Ee) .* ((sum(gNMDA.*gIntegN)).* 1./(1+exp(-0.062*gV)*Mg/3.57));
        gI = Ja + Jn;
          
        
        mIG = mIntegG ./ exp(TS./tauG_m);
        mIA = mIntegA ./ exp(TS./tauA_m);
        mIN = mIntegN + TS*(-mIntegN./tauNd_m + alphaM*mX.*(1-mIntegN));
        Xm = mX ./ exp(TS./tauNr_m);
        
        gIA = gIntegA ./ exp(TS./tauA_g);
        gIN = gIntegN + TS*(-gIntegN./tauNd_g + alphaG*gX.*(1-gIntegN));
        Xg = gX ./ exp(TS./tauNr_g); 
end


function [mitV, mitU, graV, graU] =  izhikevich(mV, mU, gV, gU, mI, gI, mFired, gFired, TS,...
                a_m, b_m, c_m, d_m, k_m, a_g, b_g, c_g, d_g, k_g, mVr, mVt, gVr, gVt, Cm, Cg)
        gV(gFired) = c_g(gFired);
        gU(gFired) = gU(gFired) + d_g(gFired);
        mV(mFired) = c_m(mFired);
        mU(mFired) = mU(mFired) + d_m(mFired);

        %Updating values for all mitral and granule cells based on
        %Izhikevich functions
        mV = mV + TS * (k_m./Cm .* (mV-mVr).*(mV-mVt) - mU./Cm + mI./Cm);
        mU = mU + TS * a_m.*(b_m .* (mV-mVr)-mU);

        gV = gV + TS * (k_g./Cg .* (gV-gVr).*(gV-gVt) - gU./Cg + gI./Cg);
        gU = gU + TS * a_g.*(b_g .* (gV-gVr)-gU);
     
        mitV = mV; mitU = mU;  
        graV = gV; graU = gU; 
end
