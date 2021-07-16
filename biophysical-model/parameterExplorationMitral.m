% Explores MC intrinsic parameters for current input
time = 1000;
TS = 0.1; tspan = 0:TS:time; 

% Setting parameters of Izhikevich model to simulate mitral cells
vr = -58;
vt = -49;

a=0.02; 
b=12; 
c=-70; 
d=13;
k = 2.5;  
Cm = 191; 

mIArray = 0:1:400;
mFreqArray = zeros(1, length(mIArray));
f = 2/1000;
Imax = 300;

delay = 0;

mVolt = zeros(1,length(tspan));

for m = 1:length(mIArray)
    spikes = 0;
    mV = vr;
    mU = 0;
    for t=1:length(tspan)
        if mV >= 30
            mV = c;
            mU = mU + d;
            if t > 0
                spikes = spikes + 1;
            end
        end
       mV = mV + TS*(k/Cm*(mV-vr)*(mV-vt) - mU/Cm + mIArray(m)/Cm);

        mU = mU + TS*a*(b*(mV-vr)-mU);
        if t > 0 && m == 200
            mVolt(t) = mV;
        end
    end
    freq = spikes / (time/1000);
    mFreqArray(1,m) = freq;
end


        