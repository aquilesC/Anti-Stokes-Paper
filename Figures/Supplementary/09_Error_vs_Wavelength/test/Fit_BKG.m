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
set(gcf,'position',[5   630   560   480])

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
    
    plot(w(:,n),s(:,n)./max(s(:,n)),'linewidth',lw)
    hold all
    ley{n} = strcat('NR',num2str(NR(n)));
end

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

%% cut the spectra at the filter
W = w(:,1);
lcut=568;
ind = find(W>lcut);

%% analize one NR to get the BKG
N = find(NR==4);
x = W(ind);
y = s(ind,N);
% defin BKG function
B = @(p,x) p(1).*exp(-(x-lcut)./p(2));
p0 = [0.2 20];
L = @(p,x) (p(1).*(p(3)/2).^2)./( (x-p(2)).^2 + (p(3)/2).^2 );
[MM,aux] = max(y);
p00= [MM x(aux) 40];
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
% plot(x,fun([p0,p00],x),'--','linewidth',lw)
xlim([550 800])
grid on

% fit
[b,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(x,y,fun_cal,[p0,p00]);

% plot fit
plot(x,fun_cal(b,x),':k','linewidth',lw)
plot(x,B([b(1),b(2)],x),'linewidth',lw)
legend('Data','Complete fit','BKG extracted')
ylabel('Norm intensity')
xlabel('Wavelength [nm]')
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
xlabel('Wavelength [nm]')

set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
% ylim([-0.05 0.05])

disp(strcat('PL BKG Extracted from NR=',num2str(NR(N)),'. b=',num2str(b)))

%% now I substract the BKG for each NR and fit only the lorentzian

% new function
fun= @(p,t) B([p(1) b(2)],x)+L([p(2) p(3) p(4)],x);

figure(4)
set(gcf,'position',[1245   634   560   480])
clf
for N=1:size(s,2)
    clear y
    y = s(ind,N);
    
    
    pl(N) = plot(x,y,'linewidth',lw);
    hold all
    
    % get the estimation for the SPR
    [MM,aux] = max(y);
    IM(N) = MM;
    [beta_aux,R2,J2,CovB2,MSE2,ErrorModelInfo2]=nlinfit(x,y,fun,[MM/10 MM*9/10 x(aux) 40]);
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
   
end

for N=1:size(s,2)
     plot(x,fun(beta(N,:),x),'--','linewidth',lw,'Color',pl(N).Color)
end

xlabel('Wavelength [nm]')
% ylim([-0.05 1.05])
% xlim([550 800])
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')


%% Show fit results
figure(5)
clf
set(gcf,'position',[3279          94         560        1022])
%
subplot(3,1,1)
errorbar(NR,SPR,eSPR,'o','MarkerSize',10,'linewidth',lw)
ylabel('SPR [nm]')
% ylim([-0.05 1.05])
% xlim([550 800])
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
%
subplot(3,1,2)
errorbar(NR,Wd,eWd,'o','MarkerSize',10,'linewidth',lw)
ylabel('SPR width [nm]')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
%
subplot(3,1,3)
errorbar(NR,A2,eA2,'o','MarkerSize',10,'linewidth',lw)
% hold all
% plot(NR,IM,'s')
ylabel('SPR amplitude [nm]')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

