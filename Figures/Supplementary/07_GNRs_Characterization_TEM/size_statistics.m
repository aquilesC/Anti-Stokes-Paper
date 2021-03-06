close all
clc
clear
%% to plot
FS =18;
BW = 2;
lw = 1.8;
%% data
d = csvread('Results.csv',1,1)
LEN =d(:,6);

for n = 1:length(LEN)/2;
    l(n) = LEN(2*n-1);
    w(n) = LEN(2*n);
end
AR = l./w;

disp(strcat('AR=',num2str(mean(AR)),'+/-',' ',num2str(std(AR))))
disp(strcat('L=',num2str(mean(l)),'+/-',' ',num2str(std(l))))
disp(strcat('W=',num2str(mean(w)),'+/-',' ',num2str(std(w))))

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
%% plot
% linecolor = [192,192,192]/255;
linecolor = [0,0,0]/255;


figure(1)
clf
set(gcf,'position',[80   790   717   308])
subplot(1,2,1)
b2 = bar(cw,nw);
set(b2,'FaceColor',[0.8500    0.3250    0.0980],'EdgeColor',linecolor,'LineWidth',lw);
hold all

xlim([10 60])
xticks(10:10:60)
xlabel('Width [nm]')
ylabel('# Nanorods')
% legend('Length','Width')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

subplot(1,2,2)
b1 = bar(cl,nl);
set(b1,'FaceColor',[0    0.4470    0.7410],'EdgeColor',linecolor,'LineWidth',lw);

xlim([20 70])
xticks(20:10:70)
xlabel('Length [nm]')
ylabel('# Nanorods')
% legend('Length','Width')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')


% figure(2)
% clf
% bar(car,nar,'k')