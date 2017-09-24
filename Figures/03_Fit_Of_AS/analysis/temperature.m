close all
clear
clc

a = load('T_01.txt');
I1 = a(:,1);
T1 = a(:,2);

a = load('T_02.txt');
I2 = a(:,1);
T2 = a(:,2);

a = load('T_03.txt');
I3 = a(:,1);
T3 = a(:,2);

% see the data
figure(1)
plot(I1,T1,'sr')
hold all
grid on
plot(I2,T2,'x','color',[0.1 0.5 0.9])
plot(I3,T3,'or')
ylim([380 500])

%%
t1 = T1-mean(T1);
t2 = T2-mean(T2);
t3 = T3-mean(T3);
% 
% figure
% plot(I1,t1,'o')
% hold all
% plot(I2,t2,'s')
% plot(I3,t3,'x')
%% fit
I = vertcat(I1,I2,I3);
T = vertcat(t1,t2,t3);
[I,ind] = sort(I)
T = T(ind);
clear ind

[p,S] = polyfit(I,T,1)
Rinv = inv(S.R)
C = (Rinv*Rinv')*S.normr^2/S.df

[y,delta] = polyval(p,0,S)

figure(2)
clf
plot(I,T,'s')
hold all
plot(I,polyval(p,I))

%% with this slope get the best origin
p1 = [p(1) mean(T1-p(1).*I1)]; % best fit with fixed slope p(1) for T1
p2 = [p(1) mean(T2-p(1).*I2)]; % best fit with fixed slope p(1) for T2
p3 = [p(1) mean(T3-p(1).*I3)]; % best fit with fixed slope p(1) for T3


[y1 delta1]=polyval(p1,0,S)
[y2 delta2]=polyval(p2,0,S)
[y3 delta3]=polyval(p3,0,S)


% %%
% x = -max(I):1:max(I);
% [yaux d] = polyval(p1,x,S);
% figure
% plot(x,yaux)
% hold all
% plot(x,yaux+d)
% plot(x,yaux-d)

x = 100:10:250;
ms = 10;
% plot
figure(3)
plot(I1,T1,'or','color',[255 165 0]/255,'linewidth',2,'markersize',ms)
hold all
plot(I2,T2,'x','color',[0 204 0]/255,'linewidth',2,'markersize',ms)
plot(I3,T3,'sr','color',[0 0 204]/255,'linewidth',2,'markersize',ms)
ylim([380 500])
plot(x,polyval(p1,x),'--','color',[255 165 0]/255,'linewidth',2)
plot(x,polyval(p2,x),'--','color',[0 204 0]/255,'linewidth',2)
plot(x,polyval(p3,x),'--','color',[0 0 204]/255,'linewidth',2)
xlabel('Laser power [\muW]','FontSize',20)
ylabel('Temperature [T]','FontSize',20)
set(gca,'linewidth',3)
set(gca,'fontsize',16)
lll=legend('20 ^oC','40 ^oC','60 ^oC');
set(lll,'location','northwest')

%% 
Treal = [293 313 333];
Tmea = [p1(2) p2(2) p3(2)]
sTmea = 7*ones(size(Tmea));

v = [283 343];
figure(4)
clf
errorbar(Treal,Tmea,sTmea,'o','linewidth',2,'markersize',ms)
hold all
plot(v,v,'linewidth',2)
xlim(v)
ylim([280 360])
xlabel('Temperature flow cell [T]','FontSize',20)
ylabel('Extracted Temperature [T]','FontSize',20)
set(gca,'linewidth',3)
set(gca,'fontsize',16)


