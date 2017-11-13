close all
clc
clear

l = linspace(600,635,10000)*1e-9; % wavelenght in m
laser = 633e-9; % m
h = 6.626e-34; % J s
kb = 1.3806e-23; %J/K
c = 3e8; % m/s
T = 273+linspace(20,1000,2);% K

E = h*c./l-h*c./laser; % J

for N =1:size(T,2)
n(:,N) = (exp(E./(kb*T(N)))-1).^-1;
n2(:,N) = (exp(E./(kb*T(N)))-1).^-1+1;
end

figure
loglog(l*1e9,n)
hold all
plot(l*1e9,n2,'--')
xlim([600 635])
ylim([1e-1 1e5])
plot([laser laser]*1e9,[1e-2 1e5],'-k')
grid on
xlabel('Wavelength [nm]')