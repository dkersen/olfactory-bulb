% Explores MC intrinsic parameters for Poisson input
time = 1000;
TS = 0.1; tspan = 0:TS:time; 


vr = -58;
vt = -49;

a=0.02; 
b=12; 
c=-70; 
d=13;
k = 2.5;  
Cm = 191;
trialNum = 50;

rateArray = 0:0.1:4;
mFreqArray = zeros(trialNum, length(rateArray));

gA = 6.7;
gN = 12;
Mg = 1;
tauA = 14.3108;
tauNr = 13;
tauNd = 70;
f = 6/1000;
alpha = 0.03;

numRecN = 100;
numRecA = 100;
mIAHist = zeros(1,length(tspan));
mINHist = zeros(1,length(tspan));
mIHist = zeros(1, length(tspan));
W = 0.5;

for r = 1:length(rateArray)
    for trial = 1:trialNum
        spikes = 0;
         mV = vr;
        mU = 0;
        mIntegN = zeros(1, numRecN);
        mIntegA = zeros(1, numRecA);
        mX = zeros(1, numRecN);
        Ee = 0;
        for t=1:length(tspan)
            if mV >= 30
                mV = c;
                mU = mU + d;
                spikes = spikes + 1;
            end
            R = (rand(size(mIntegA)) < rateArray(r)/3 * TS/1000 + rateArray(r)/3 * TS/1000 * (sin(2*pi*f*tspan(t) - pi/2) + 1));
            mIntegA = mIntegA + R*W.*(1-mIntegA);
            mX = mX + R*W.*(1-mX);
            
            mI = -gA * sum(mIntegA) * (mV-Ee) - gN * sum(mIntegN)*(mV - Ee) .* 1./(1+exp(-0.062*mV)*Mg/3.57);
%             mIAHist(t) = mIntegA; %-gA * sum(mIntegA) * (mV-Ee);
%             mINHist(t) = mIntegN; %-gN * sum(mIntegN)*(mV - Ee) .* 1./(1+exp(-0.062*mV)*Mg/3.57);
%             mIHist(t) = mI;
            mV = mV + TS*(k/Cm*(mV-vr)*(mV-vt) - mU/Cm + mI/Cm);
            mU = mU + TS*a*(b*(mV-vr)-mU);
            mIntegA = mIntegA/exp(TS/tauA);
            mIntegN = mIntegN + TS*(-mIntegN./tauNd + alpha*mX.*(1-mIntegN));
            mX = mX/exp(TS/tauNr);
        end
        mFreqArray(trial,r) = spikes / (time/1000);
    end
end


        