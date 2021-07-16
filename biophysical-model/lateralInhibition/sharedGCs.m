% Plot the mean number of shared GCs between MCs as a function of
% distance, as well as the fit

load('misc.mat')
load('mitPairs.mat')

x = zeros(1,length(mitPairs)); 
y = zeros(1,length(mitPairs));
for i = 1:length(mitPairs)
    x(i) = distMat(mitPairs(i,1), mitPairs(i,2));
    y(i) = overlapMat(mitPairs(i,1), mitPairs(i,2));
end

interval=100;
edges = 0:interval:ceil(max(x)/interval)*interval;

meanz = zeros(1,length(edges)-1);
errz = zeros(1, length(edges)-1);

for i = 1:length(edges)-1
    Z = find(x >= edges(i) & x < edges(i+1));
    meanz(i) = mean(y(Z));
    errz(i) = std(y(Z))/sqrt(length(Z));
end

xmid = 0.5*(edges(1:end-1)+edges(2:end));

figure
hold on

x = 0:1:1200;
a1 = 229.2;
b = -0.0001721;
n = 1.545;
fit = a1*exp(b*x.^n); 
%r^2 = 0.9997
% a1 =       212.1;
% c1 =       290.2;
% fit = a1.*exp(-((x)./c1).^2);
errorbar(xmid, meanz, errz,'o')
plot(x,fit)
