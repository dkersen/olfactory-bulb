% Explores GC intrinsic parameters for current input

time = 1000;
TS = 0.1; tspan = 0:TS:time; 

% Setting parameters of Izhikevich model to simulate granule cells
vr = -71;
vt = -39;

a=0.01; 
b=-2/15; 
c=-75; 
d=1.2;
k = 1/15;  
Cg = 48;

gIArray = 0:1:100;
nIArray = 0;
gFreqArray = zeros(length(gIArray),length(nIArray));
f=6/1000;

delay = 0;

gVolt = zeros(1,length(tspan));

for g = 1:length(gIArray)
    for n = 1:length(nIArray)
        spikes = 0;
        gV = vr;
        gU = 0;
        for t=1:length(tspan)
            if gV >= 25
                gV = c;
                gU = gU + d;
                if t > 0
                    spikes = spikes + 1;
                end
            end
            gV = gV + TS*(k/Cg*(gV-vr)*(gV-vt) - gU/Cg + gIArray(g)/Cg);
            gU = gU + TS*a*(b*(gV-vr)-gU);
            if t > 0 && gIArray(g) == 45
                gVolt(t) = gV;
            end
        end
        freq = spikes / (time/1000);
        gFreqArray(g,n) = freq;
    end
end

        
