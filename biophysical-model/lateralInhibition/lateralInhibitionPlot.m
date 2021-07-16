% Plot the strength of lateral inhibition as a function of distance and fit

load('differences')
load('distanceVector')

edges = 0:100:1200;

x = distVec;
y = diff;

meanz = zeros(1,length(edges)-1);
errz = zeros(1, length(edges)-1);

for i = 1:length(edges)-1
    Z = find(x >= edges(i) & x < edges(i+1));
    meanz(i) = mean(y(Z));
    errz(i) = std(y(Z))/sqrt(length(Z));
end

xmid = 0.5*(edges(1:end-1)+edges(2:end));

figure
x = 0:1:1200;
a1 = 11.08;
b = -9.375e-06;
n= 1.976; %0.9972
fit = a1*exp(b*x.^n);

errorbar(xmid, meanz, errz,'o')
hold on
plot(x,fit);



