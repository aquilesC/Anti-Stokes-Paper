close all
clear
clc

%% to plot
FS =18;
BW = 2;
lw = 1.8;

%% load data

f = dir('532nm_spectra*.txt');

figure(1)
set(gcf,'position',[5   400   560   700])
subplot(2,1,1)
for n=1:length(f)
    a = load(f(n).name);
    
    aux = strfind(f(n).name,'NR');
    NR(n) =str2num(f(n).name(aux(1)+2:end-4));
    clear aux
    %     if n==1
    %         w(:,n) = a(1,:)';
    %         s(:,n) = a(2,:)';
    %     else
    %         w(:,n) = a(:,1);
    %         s(:,n) = a(:,2);
    %     end
    w(:,n) = a(:,1);
    s(:,n) = a(:,2);
    clear a
    figure(1)
    subplot(2,1,1)
    plot(w(:,n),s(:,n)./max(s(:,n)),'linewidth',lw)
    hold all
    
    subplot(2,1,2)
    plot(1239.8./w(:,n),s(:,n)./max(s(:,n)),'linewidth',lw)
    hold all
    ley{n} = strcat('NR =',num2str(NR(n)));
end
subplot(2,1,1)
ylim([-0.05 1.05])
xlim([550 800])
legend(ley)
ylabel('Norm intensity')
xlabel('Wavelength [nm]')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')


subplot(2,1,2)
ylim([-0.05 1.05])
xlim([1.5 2.5])
legend(ley)
ylabel('Norm intensity')
xlabel('E [eV]')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

%% cut the spectra at the filter
W = w(:,1); % wavelegnth in nm
E = 1239.8./W; % energy in eV
lcut=568;
ind = find(W>lcut);

%% analize one NR to get the BKG
N = find(NR==4);
x = E(ind);
y = s(ind,N);
% defin BKG function
B = @(p,x) p(1).*exp(x./p(2));
p0 = [1e-5 0.1];
L = @(p,x) (p(1).*(p(3)/2).^2)./( (x-p(2)).^2 + (p(3)/2).^2 );
[MM,aux] = max(y);
p00= [MM x(aux) 0.15];
clear aux MM

fun_cal = @(p,t) B([p(1) p(2)],x)+L([p(3) p(4) p(5)],x);

% plot initial guess
figure(2)
set(gcf,'position',[680   634   560   480])
clf
plot(x,y,'linewidth',lw)
hold all
% plot(x,B(p0,x),'linewidth',lw)
% plot(x,L(p00,x),'linewidth',lw)
% plot(x,fun_cal([p0,p00],x),'--','linewidth',lw)

grid on

% fit
[b,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(x,y,fun_cal,[p0,p00]);

% plot fit
plot(x,fun_cal(b,x),':k','linewidth',lw)
plot(x,B([b(1),b(2)],x),'linewidth',lw)
legend('Data','Complete fit','BKG extracted','Location','NorthWest')
ylabel('Norm intensity')
xlabel('E [eV]')
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

% plot res


figure(3)
clf
set(gcf,'position',[680   290   560   250])
plot(x,R,'linewidth',lw)
ylabel('Residuals')
grid on
xlabel('E [eV]')

set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
% ylim([-0.05 0.05])

disp(strcat('PL BKG Extracted from NR=',num2str(NR(N)),'. b=',num2str(b)))
% with this for all the NR I create the file BKG_E_results.txt
%%
a = load('BKG_E_results.txt');

D = mean(a(:,3));
eD = std(a(:,3));
clear a

%% now get the SPR
% new function
fun   = @(p,x) B([p(1) D],x)+L([p(2) p(3) p(4)],x);
fun_u = @(p,x) B([p(1) D+eD/2],x)+L([p(2) p(3) p(4)],x);
fun_d = @(p,x) B([p(1) D-eD/2],x)+L([p(2) p(3) p(4)],x);

% figure(4)
% set(gcf,'position',[1245   634   560   480])
% clf
for N=1:size(s,2)
    clear y
    y = s(ind,N);
    figure(100+N)
    clf
    pl(N) = plot(x,y,'linewidth',lw);
    hold all
    disp('Center')
    % get the estimation for the SPR
    [MM,aux] = max(y);
    IM(N) = MM;
    [beta_aux,R2,J2,CovB2,MSE2,ErrorModelInfo2]=nlinfit(x,y,fun,[b(1) MM x(aux) 0.15]);
    clear MM aux
    
    beta(N,:) = beta_aux;
    Res(N,:) = R2;
    Co(N,:,:) = CovB2;
    
    A1(N) = beta_aux(1);
    A2(N) = beta_aux(2);
    SPR(N) = beta_aux(3);
    Wd(N) = beta_aux(4);
    
    eA1(N) = sqrt(CovB2(1,1));
    eA2(N) = sqrt(CovB2(2,2));
    eSPR(N) = sqrt(CovB2(3,3));
    eWd(N) = sqrt(CovB2(4,4));
    
    clear beta_aux
    
    disp('UP')
    % get the estimation for the SPR upper bound
    [MM,aux] = max(y);
    IM(N) = MM;
    [beta_aux,R2,J2,CovB2,MSE2,ErrorModelInfo2]=nlinfit(x,y,fun_u,[b(1) MM x(aux) 0.15]);
    clear MM aux
    
    beta_u(N,:) = beta_aux;
    clear beta_aux
    
    disp('DOWN')
    % get the estimation for the SPR upper bound
    [MM,aux] = max(y);
    IM(N) = MM;
    [beta_aux,R2,J2,CovB2,MSE2,ErrorModelInfo2]=nlinfit(x,y,fun_d,[b(1) MM x(aux) 0.15]);
    clear MM aux
    
    beta_d(N,:) = beta_aux;
    clear beta_aux
    
    
end

for N=1:size(s,2)
    figure(100+N)
    plot(x,fun(beta(N,:),x),'--r','linewidth',lw)
    plot(x,fun_u(beta_u(N,:),x),':r','linewidth',lw-1)
    plot(x,fun_d(beta_d(N,:),x),':r','linewidth',lw-1)
    
    text(x(end),max(fun(beta(N,:),x))*0.9,strcat('SPR=',num2str(SPR(n))))
    text(x(end),max(fun(beta(N,:),x))*0.7,strcat('SPR Width=',num2str(Wd(n))))
    ylim([0 max(fun(beta(N,:),x))*1.1])
    title(strcat('NR=',num2str(NR(N))))
    xlabel('E [eV]')
    % ylim([-0.05 1.05])
    % xlim([550 800])
    grid on
    set(gca,'FontSize',FS)
    set(gca,'Linewidth',BW)
    set(gca,'XMinorGrid','off')
    set(gca,'YMinorGrid','off')
    
end

%% prepare figure for SI
to_plot = [4 78 33];

figure(4)
clf
set(gcf,'position',[408   329   647   487])

for n=1:length(to_plot)
    clear a
    figure(4)
    plot(E,s(:,find(NR==to_plot(n)))./max(s(:,find(NR==to_plot(n)))),'linewidth',lw)
    hold all
    
end


xlim([1.6 2.4])
ylim([0 1])
plot([1239.8/532 1239.8/532],[0 1],'--','linewidth',lw,'color',[0,153,51]/255)
% ylabel('Norm intensity')
% xlabel('Wavelength [nm]')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

% saveas(gcf,'Normalized_spectra_532nm.fig','fig')
% saveas(gcf,'Normalized_spectra_532nm.pdf','pdf')

