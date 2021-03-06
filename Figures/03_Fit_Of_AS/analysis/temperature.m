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

%% fit
I = vertcat(I1,I2,I3);
T = vertcat(t1,t2,t3);
[I,ind] = sort(I)
T = T(ind);
clear ind

[p,S] = polyfit(I,T,1);
Rinv = inv(S.R);
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
clf
set(gcf,'units','centimeters')
set(gcf,'position',[10 20 85*3 55*3]/10)
% plot the fits
plot(x,polyval(p1,x),'--','color',[255 165 0]/255,'linewidth',2)
hold all
plot(x,polyval(p2,x),'--','color',[0 204 0]/255,'linewidth',2)
plot(x,polyval(p3,x),'--','color',[0 0 204]/255,'linewidth',2)

errores = true;
if errores== false
    hh1 = plot(I1,T1,'or','color',[255 165 0]/255,'linewidth',2,'markersize',ms,'markerfacecolor','w');
    hh2 = plot(I2,T2,'^','color',[0 204 0]/255,'linewidth',2,'markersize',ms,'markerfacecolor','w');
    hh3 = plot(I3,T3,'sr','color',[0 0 204]/255,'linewidth',2,'markersize',ms,'markerfacecolor','w');
else
    hh1 = errorbar(I1,T1,ones(size(I1)).*5,'or','color',[255 165 0]/255,'linewidth',2,'markersize',ms,'markerfacecolor','w');
    hh2 = errorbar(I2,T2,ones(size(I2)).*5,'^','color',[0 204 0]/255,'linewidth',2,'markersize',ms,'markerfacecolor','w');
    hh3 = errorbar(I3,T3,ones(size(I3)).*5,'sr','color',[0 0 204]/255,'linewidth',2,'markersize',ms,'markerfacecolor','w');
end
ylim([380 500])

% xlabel('Laser power [\muW]','FontSize',20)
% ylabel('Temperature [T]','FontSize',20)
set(gca,'linewidth',3)
set(gca,'fontsize',16)
lll=legend([hh1 hh2 hh3],{'20 ^oC','40 ^oC','60 ^oC'});
set(lll,'location','northwest')

% saveas(gcf,'TvsInt','fig')
% saveas(gcf,'TvsInt','pdf')


%% 
Treal = [293 313 333];
Tmea = [p1(2) p2(2) p3(2)];
sTmea = 8*ones(size(Tmea));



v = [280 343];
figure(4)
set(gcf,'units','centimeters')
set(gcf,'position',[10 20 25*5 18*5]/10)
clf
errorbar(Treal,Tmea,sTmea/2,'s','linewidth',2,'color',[0 0 204]/255,'markersize',ms,'MarkerEdgeColor',[0 0 204]/255,...
    'MarkerFaceColor',[0 0 204]/255)
hold all
plot(v,v,'r','linewidth',2)
xlim(v)
ylim([280 360])
plot(v,v+ori,'--','linewidth',2')
% xlabel('Temperature flow cell [T]','FontSize',20)
% ylabel('Extracted Temperature [T]','FontSize',20)
set(gca,'linewidth',3)
set(gca,'fontsize',16)
xticks(280:20:340)
yticks(280:20:360)

% saveas(gcf,'Inset','fig')
% saveas(gcf,'Inset','pdf')
%% from aquiles
p1a =[0.71225756  298.08501504];
p2a =[0.71225756  316.62098401];
p3a =[0.71225756  340.344689644];
% Tmea_a = [298.08501504  316.62098401  340.34468964]; %witout specifiy tol
% with tol 1e-3 in the jupyter program
Tmea_a = [296  315  339]
sTmea_a = 6*ones(size(Tmea_a));
ori = mean(Tmea_a-Treal);

Treal

figure(5)
set(gcf,'units','centimeters')
set(gcf,'position',[10 20 25*5 18*5]/10)
clf
plot(v,v+ori,'--','linewidth',2','color',[0 0 204]/255)
hold all
plot(v,v,'r','linewidth',2)
errorbar(Treal,Tmea_a,sTmea_a/2,'s','linewidth',2,'color',[0 0 204]/255,'markersize',ms,'MarkerEdgeColor',[0 0 204]/255,...
    'MarkerFaceColor','w')



% xlabel('Temperature flow cell [T]','FontSize',20)
% ylabel('Extracted Temperature [T]','FontSize',20)
set(gca,'linewidth',3)
set(gca,'fontsize',16)
xticks(280:20:340)
yticks(280:20:360)
xlim(v)
ylim([280 360])
saveas(gcf,'Inset_a','fig')
saveas(gcf,'Inset_a','pdf')

