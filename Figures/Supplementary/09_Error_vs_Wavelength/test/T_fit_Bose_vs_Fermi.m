clear 
close all
clc

El = 1239.8/633; % laser energy
Elim = 1239.8/530;
Kb = 8.6173303*1e-5; % eV /K
n_oc = @(p,x) p(1)./(exp((x-El)./(Kb.*p(2)))-1); % ocupation number
n_oc_f = @(p,x) p(1)./(exp((x-El)./(Kb.*p(2)))+1); % ocupation number

E = linspace(El,Elim,1000); % eV


T=300:50:500; % Temp in K

for n=1:length(T)
S_bose(n,:) = n_oc([1,T(n)],E);
S_fer(n,:) = n_oc_f([1,T(n)],E);
end

figure
semilogy(E,S_bose)
hold all
semilogy(E,S_fer,'--')

