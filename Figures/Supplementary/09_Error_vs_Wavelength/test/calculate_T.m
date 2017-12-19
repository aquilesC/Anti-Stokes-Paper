close all
clear
clc

%% to plot
FS =18;
BW = 2;
lw = 1.8;
%% load BKG info
a = load('BKG_E_results.txt');
b=mean(a(:,2));
D = mean(a(:,3));
eD = std(a(:,3));
clear a

% MOdel functions for the PL at 532nm

% BKG (E)
B = @(p,x) p(1).*exp(x./p(2));

% Plasmon: Lorentzian (E)
L = @(p,x) (p(1).*(p(3)/2).^2)./( (x-p(2)).^2 + (p(3)/2).^2 );

fun   = @(p,x) B([p(1) D],x)+L([p(2) p(3) p(4)],x);
fun_u = @(p,x) B([p(1) D+eD],x)+L([p(2) p(3) p(4)],x);
fun_d = @(p,x) B([p(1) D-eD],x)+L([p(2) p(3) p(4)],x);
%% load data

s532 = load('aS532nm_spectra_NR1.txt')';
% reduce to remove the filter
ind = find(s532(:,1)>543);
W1 = s532(ind,1);
S1 = s532(ind,2);
E1 = 1239.8./W1;
clear s532 ind


%% plot
figure(1)
clf
plot(E1,S1)
hold all
% fit SPR
[MM,aux] = max(S1);
beta0 = [b/10 MM E1(aux) 0.15];
% plot(E,fun(beta0,E))
% fit
[beta,R2,J2,CovB2,MSE2,ErrorModelInfo2]=nlinfit(E1,S1,fun,beta0);
[beta_u,R2,J2,CovB2,MSE2,ErrorModelInfo2]=nlinfit(E1,S1,fun_u,beta0);
[beta_d,R2,J2,CovB2,MSE2,ErrorModelInfo2]=nlinfit(E1,S1,fun_d,beta0);


plot(E1,fun(beta,E1),'--r','linewidth',lw)
plot(E1,fun_u(beta_u,E1),':r','linewidth',lw-1)
plot(E1,fun_d(beta_d,E1),':r','linewidth',lw-1)
xlabel('E [eV]')
ylim([0 max(S1)*1.1])
xlim([1.5 2.3])
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

text(E1(end),max(fun(beta,E1))*0.9,strcat('SPR=',num2str(beta(3)),' eV'),'FontSize',16)
text(E1(end),max(fun(beta,E1))*0.7,strcat('SPR Width=',num2str(beta(4)),' eV'),'FontSize',16)

%% load 633 data
f = dir('aS633nm_spectra*.txt');

figure(2)
set(gcf,'position',[12   230   650   800])
clf
for n=1:length(f)
    a = load(f(n).name)';
    
    aux = strfind(f(n).name,'pw');
    pw(n) =str2num(f(n).name(aux(1)+3:end-4));
    clear aux
    
    w(:,n) = a(:,1);
    s(:,n) = a(:,2);
    clear a
    figure(2)
    subplot(2,1,1)
    plot(w(:,n),s(:,n)./max(s(:,n)),'linewidth',lw)
    hold all
    subplot(2,1,2)
    plot(1239.8./w(:,n),s(:,n)./max(s(:,n)),'linewidth',lw)
    hold all
    
    ley{n} = strcat('pw =',num2str(pw(n)));
    
end
subplot(2,1,1)
ylim([-0.05 1.05])
xlim([550 650])
legend(ley,'Location','NorthWest')
ylabel('Norm intensity')
xlabel('Wavelength [nm]')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

subplot(2,1,2)
ylim([-0.05 1.05])
xlim([1239.8./650 1239.8./550])
ylabel('Norm intensity')
xlabel('Wavelength [nm]')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')



%% fit the temp
W = w(:,1);
E =1239.8./W;
% range = find(E>2.0325 & E<2.1014);
range = find(E>2.015 & E<2.102);
range = find(E>2.015 & E<2.12);


El = 1239.8/633; % laser energy
Kb = 8.6173303*1e-5; % eV /K
n_oc = @(p,x) p(1)./(exp((x-El)./(Kb.*p(2)))-1); % ocupation number

%% for to fit every spectra
for n=1:size(pw,2)
    x=E(range);
    y = s(range,n);
    ey = ones(size(y));
    %          ey = real(sqrt((y)))+1e-4;
    
    ini = [1 300];
    disp('Center')
    model = @(p,x) L([beta(2),beta(3),beta(4)],x).*n_oc([p(1) p(2)],x);
    [Beta_T,R3,J3,CovB3,MSE3,ErrorModelInfo3]=nlinfit(x,y,model,ini,'Weights',ey);
    N0(n) = Beta_T(1);
    T(n) = Beta_T(2);
    eT(n) = sqrt(CovB3(2,2));
    
%     disp('UP')
%     model = @(p,x) L([beta_u(2),beta_u(3),beta_u(4)],x).*n_oc([p(1) p(2)],x);
%     
%     [Beta_Tu,R3,J3,CovB3,MSE3,ErrorModelInfo3]=nlinfit(x,y,model,ini,'Weights',ey);
%     N0u(n) = Beta_T(1);
%     Tu(n) = Beta_T(2);
%     eTu(n) = sqrt(CovB3(2,2));
%     
%     
%     disp('Down')
%     model = @(p,x) L([beta_d(2),beta_d(3),beta_d(4)],x).*n_oc([p(1) p(2)],x);
%     
%     [Beta_Td,R3,J3,CovB3,MSE3,ErrorModelInfo3]=nlinfit(x,y,model,ini,'Weights',ey);
%     N0d(n) = Beta_T(1);
%     Td(n) = Beta_T(2);
%     eTd(n) = sqrt(CovB3(2,2));
    
    
    figure(100+n)
    clf
%     errorbar(x,y,ey,'o','MarkerSize',3)
    plot(x,y,'linewidth',lw)
    set(gca,'yscale','log')
    hold all
    plot(x,model(Beta_T,x),'--','linewidth',lw)
    
    %     plot(E,L([A2(n),SPR(n),Wd(n)],E))
    %     plot(E,s(:,n))
    %     plot([El El],[0 max(s2(:,n))],'--k')
    grid on
    %     xlim([2 2.15])
    %     ylim([1e0 1e3])
    title(strcat('PW=',num2str(pw(n)),'mW. T=',num2str(Beta_T(2))))
    ylabel('Intensity')
    xlabel('E [eV]')
    grid on
    set(gca,'FontSize',FS)
    set(gca,'Linewidth',BW)
    set(gca,'XMinorGrid','off')
    set(gca,'YMinorGrid','off')
    
    
    
    figure(200+n)
    clf
    plot(E,s(:,n),'linewidth',lw)
    hold all
    plot(x,model(Beta_T,x),'--','linewidth',lw)
    plot(E,L([beta(2),beta(3),beta(4)],E),'--','linewidth',lw)
    
    %     plot(E,L([A2(n),SPR(n),Wd(n)],E))
    %     plot(E,s(:,n))
    %     plot([El El],[0 max(s2(:,n))],'--k')
    grid on
    %     xlim([2 2.15])
    %     ylim([1e0 1e3])
    title(strcat('PW=',num2str(pw(n)),'mW. T=',num2str(Beta_T(2))))
    ylabel('Intensity')
    xlabel('E [eV]')
    grid on
    set(gca,'FontSize',FS)
    set(gca,'Linewidth',BW)
    set(gca,'XMinorGrid','off')
    set(gca,'YMinorGrid','off')
    xlim([1.9 2.5])
end


%%
% P = polyfitweighted(pw,T,1,eT);
[P,S] = polyfit(pw,T,1);
[~,delta] = polyval(P,0,S);
figure(4)
clf
errorbar(pw,T,eT,'o','MarkerSize',10,'linewidth',lw)
hold all
% errorbar(pw,Tu,eTu,'o','MarkerSize',10,'linewidth',lw)
% errorbar(pw,Td,eTd,'o','MarkerSize',10,'linewidth',lw)
plot(pw,polyval(P,pw),'-')
% plot([0 max(pw)],ones(1,2).*(23+273))
grid on
xlabel('Power [mW]')
ylabel('Temperature [K]')
title(strcat('T_{room}=',num2str(floor(P(2))),'(',num2str(floor(delta)),' )K'))
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot for the SI
close all

figure(1)
clf
set(gcf,'position',[30  400  260  260])
plot(E1,S1,'o','MarkerSize',5,'MarkerFaceColor','w')
hold all
plot(E1,fun(beta,E1),'--k','linewidth',lw)
xlabel('E [eV]')

ylim([0 max(S1)])
xlim([1.6 2.3])
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')


saveas(gcf,'PL_spectra.fig','fig')
saveas(gcf,'PL_spectra.pdf','pdf')

%%
figure(2)
clf
MAR=['^','o','s']
aux = [1,3,6];
range = find(E>2.01);
x=E(range);
for n=1:size(aux,2)
    figure(2)
    pl(n) = plot(E,s(:,aux(n)),MAR(n),'MarkerSize',6,'MarkerFaceColor','w');
%     set(gca,'yscale','log')
    hold all
    legd{n} = strcat('P = ',num2str(floor(pw(aux(n)))),' mW');
end

for n=1:size(aux,2)
       plot(x,model([N0(aux(n)) T(aux(n))],x),'-','linewidth',lw,'Color',pl(n).Color)
end

% figure
% for n=1:size(aux,2)
%     plot(x,model([N0(aux(n)) T(aux(n))],x),'--','linewidth',lw,'Color',pl(n).Color)
%     hold all
% end


ylim([0 600])
xlim([1.95 2.3])
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
legend(legd{1:length(aux)})

clear aux


% saveas(gcf,'AS_spectra.fig','fig')
% saveas(gcf,'AS_spectra.pdf','pdf')

%%
x_lim = [50 110];
figure(3)
set(gcf,'position',[966.6000  333.8000  560.0000  420.0000])
clf
plot(x_lim,polyval(P,x_lim),'-k','linewidth',lw)
hold all
errorbar(pw,T,eT,'s','MarkerSize',10,'linewidth',lw,'MarkerFaceColor','w')
set(gca,'XTick',[50:15:115])

xlim(x_lim)
grid on
xlabel('Power [mW]')
ylabel('Temperature [K]')
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

disp(strcat('T_{room}=',num2str(floor(P(2))),'(',num2str(floor(delta)),' )K'))

% saveas(gcf,'TvsPower.fig','fig')
% saveas(gcf,'TvsPower.pdf','pdf')
