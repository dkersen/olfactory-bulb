% Calculate the average LFP across trials, either for the reduced (smol)
% network or the full network

TS = 0.1;
fs = 1/(TS/1000);
fcut = 200;
[b,a] = butter(6,fcut/(fs/2),'low');

freqlength = 100; 
window_size = 4000;
shift = round(window_size/2);
numTrials = 10;

for trial = 1:numTrials
   fname = sprintf('LFPresultsSM6smol_results_trial3%d.mat',trial);
%    fname = sprintf('LFPresultsSM6_results_trial3%d.mat',trial);
    load(fname)
    
    LFP_tot = LFP_NMDA+LFP_AMPA+LFP_GABA;
    y = filtfilt(b,a,LFP_tot);

    y1 = y(2001:10000);
    y1 = detrend(y1);
    y2 = y(12001:end);
    y2 = detrend(y2);


    [pyy, f] = pwelch(y1,window_size,shift,[],fs);

    [pyy2, f] = pwelch(y2,window_size,shift,[],fs);
    
    if trial == 1
        pyy_tot = zeros(numTrials,length(pyy));
        pyy_tot_2 = zeros(numTrials,length(pyy2));
    end
    
    pyy_tot(trial,:) = pyy;
    pyy_tot_2(trial,:) = pyy2;
end

f = f';
LFP_mean = mean(pyy_tot);
LFP_SEM = std(pyy_tot)/sqrt(numTrials);
LFP_mean_2 = mean(pyy_tot_2);
LFP_SEM_2 = std(pyy_tot_2)/sqrt(numTrials);

% errorbar(f(1:freqlength),LFP_mean(1:freqlength),LFP_SEM(1:freqlength))
% hold on
% errorbar(f(1:freqlength),LFP_mean_2(1:freqlength),LFP_SEM_2(1:freqlength))
% xlim([0,100])

figure(1)
errorbar(f(1:freqlength),LFP_mean(1:freqlength),LFP_SEM(1:freqlength))
xlim([0,100])

figure(2)
errorbar(f(1:freqlength),LFP_mean_2(1:freqlength),LFP_SEM_2(1:freqlength))
xlim([0,100])