close all
clc
clear

d = csvread('Results.csv',1,1)
LEN =d(:,6);

for n = 1:length(LEN)/2;
    l(n) = LEN(2*n-1);
    w(n) = LEN(2*n);
end
AR = l./w;

median(AR)
median(l)
median(w)
%%
NBIN = 6;
[a1 a2] = histcounts(l,NBIN );
nl = a1;
cl = a2(1)+mean(diff(a2))/2:mean(diff(a2)):max(a2);
clear a1 a2

[a1 a2] = histcounts(w,NBIN );
nw = a1;
cw = a2(1)+mean(diff(a2))/2:mean(diff(a2)):max(a2);
clear a1 a2


[a1 a2] = histcounts(AR,NBIN );
nar = a1;
car = a2(1)+mean(diff(a2))/2:mean(diff(a2)):max(a2);
clear a1 a2

figure(1)
clf
bar(cl,nl)
hold all
bar(cw,nw,'r')

figure(2)
clf
bar(car,nar,'k')