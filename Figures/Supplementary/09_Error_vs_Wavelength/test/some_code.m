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
   
end

for N=1:size(s,2)
     plot(x,fun(beta(N,:),x),'--','linewidth',lw,'Color',pl(N).Color)
end

xlabel('E [eV]')
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
ylabel('SPR [eV]')
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
ylabel('SPR width [eV]')
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
ylabel('SPR amplitude')
grid on
set(gca,'FontSize',FS)
set(gca,'Linewidth',BW)
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
%% use the 633 data to fit the temp

% load 633 data
f2 = dir('633nm_spectra*.txt');

figure(6)
set(gcf,'position',[500   400   560   700])
subplot(2,1,1)
for n=1:length(f2)
    a = load(f2(n).name);
    
    aux = strfind(f2(n).name,'NR');
    NR2(n) =str2num(f2(n).name(aux(1)+2:end-4));
    clear aux
    %     if n==1
    %         w(:,n) = a(1,:)';
    %         s(:,n) = a(2,:)';
    %     else
    %         w(:,n) = a(:,1);
    %         s(:,n) = a(:,2);
    %     end
    w2(:,n) = a(:,1);
    s2(:,n) = a(:,2);
    clear a
    figure(6)
    subplot(2,1,1)
    plot(w2(:,n),s2(:,n)./max(s2(:,n)),'linewidth',lw)
    hold all
    
    subplot(2,1,2)
    plot(1239.8./w2(:,n),s2(:,n)./max(s2(:,n)),'linewidth',lw)
    hold all
    ley{n} = strcat('NR =',num2str(NR2(n)));
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

%% select the range and model
close all

range = find(E>2.01 & E<2.2);

El = 1239.8/633; % laser energy
Kb = 8.6173303*1e-5; % eV /K 
n_oc = @(p,x) p(1)./(exp((x-El)./(Kb.*p(2)))-1); % ocupation number


% figure(7)
% clf
% set(gcf,'position',[900   400   560   700])
for n=1:size(NR,2)
    x=E(range);
    y = s2(range,n);
    ey = sqrt(E(range));
    
    ini = [1 300];
    model = @(p,x) L([A2(n),SPR(n),Wd(n)],x).*n_oc([p(1) p(2)],x);
    
    [Beta_T,R3,J3,CovB3,MSE3,ErrorModelInfo3]=nlinfit(x,y,model,ini);
    T0(n) = Beta_T(1);
    T(n) = Beta_T(2);
    eT(n) = sqrt(CovB3(2,2)) 
    figure()
    plot(E,s2(:,n),'o','MarkerSize',3)
%     set(gca,'yscale','log')
    hold all
    plot(x,model(Beta_T,x),'--')
    
    plot(E,L([A2(n),SPR(n),Wd(n)],E))
    plot(E,s(:,n))
    plot([El El],[0 max(s2(:,n))],'--k')
    grid on
%     xlim([2 2.15])
%     ylim([1e0 1e3])
    title(strcat('NR=',num2str(NR(n)),'. Fit=',num2str(Beta_T)))
end


%% 
figure
errorbar(SPR,T,eT,'*')
hold all
plot([El El],[min(T) max(T)],'--')