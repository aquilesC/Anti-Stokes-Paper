%% from aquiles code
clear all
clc
close all
%%
data_folder = '\\data02\pi-orrit\Aquiles\Data\Background-Free\2014-06-05\'
name_532 = 'S140604_532';
name_633 = 'S140604_633_01';
%% import the coordinates from the found particles
fid = fopen([data_folder name_532 '_final.txt'],'r');
fgetl(fid);
data_532 = textscan(fid,'%*s%f%f%*f%*s','Delimiter',',');
fclose(fid);
data_532 = [data_532{1} data_532{2}];
num_bkg = 4;
data_532(end-num_bkg+1:end,:) = [];
fid = fopen([data_folder name_633 '_final.txt'],'r');
fgetl(fid);
data_633 = textscan(fid,'%*s%f%f%*f%*s%*s','Delimiter',',');
fclose(fid);
data_633 = [data_633{1} data_633{2}];
num_bkg = 3;
data_633(end-num_bkg+1:end,:) = [];
%% Iterate to find minimum distance between the 532 and 633 scans
conversion = zeros(2,size(data_633,1)); %Which particle in the 532 scan corresponds to which in the 633

for i = 1:size(data_633,1)
    coord_633 = data_633(i,:);
    for j=1:size(data_532,1)
        coord_532 = data_532(j,:);
        dist_new = sqrt(sum((coord_532-coord_633).^2));
        if j==1 || dist_new <= dist
            dist = dist_new;
            conversion(1,i) = j;
            conversion(2,i) = dist;
        end
    end
end

%% Get the 532 spectra
spectra_to_process = 'AuNR_S140604_532_1.txt';
number_of_backgrounds = 4;
% Read File (ASCII)
spectra = importdata(strcat(data_folder,spectra_to_process)); 

% Prepare data for processing
lambda = spectra(1:1340,1); % Extract wavelength information (only once)
spectra(:,1) = []; % Erase wavelength from the matrix

spectra = reshape(spectra', 1340,length(spectra)/1340);
spectra_bkg = spectra(:,end-number_of_backgrounds+1:end);
spectra(:,end-number_of_backgrounds+1:end) = [];
spectra_bkg = sum(spectra_bkg,2)/number_of_backgrounds;
spectra_532 = spectra;
spectra_bkg_532 = spectra_bkg;

%% Process spectra of the 633
spectra_to_process = 'AuNR_S140604_633_1.txt';
number_backgrounds = 3;

% Read File (ASCII)
spectra = importdata(strcat(data_folder,spectra_to_process)); 

% Prepare data for processing
lambda = spectra(1:1340,1); % Extract wavelength information (only once)
spectra(:,1) = []; % Erase wavelength from the matrix

spectra = reshape(spectra', 1340,length(spectra)/1340);
spectra_bkg = spectra(:,end-number_backgrounds+1:end);
spectra(:,end-number_backgrounds+1:end) = [];
spectra_bkg = sum(spectra_bkg,2)/number_backgrounds;
spectra1 = spectra;
background1 = spectra_bkg;


%% TO EXPORT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% find spr
SPR = zeros(1,size(spectra1,2));
figure(1)
clf
for i = 1:size(spectra1,2)
   to_fit = spectra_532(:,conversion(1,i))-spectra_bkg_532;
   S(i,:) = to_fit;
   [ind1 ind2] =max(to_fit);
   SPR(i) = lambda(ind2);   
   figure(1)
   plot(lambda,to_fit)
   hold all
end
%
figure
plot(SPR)

%% select the particle and save
i=78
spectra_export = S(i,:);

figure
plot(lambda,spectra_export)
% to save 532

A = [lambda,spectra_export'];
save_folder = 'D:\Martin\Publications\anti-stokes\Figures\Supplementary\09_Error_vs_Wavelength\test\';
fileID = fopen(strcat(save_folder,'532nm_spectra_NR',num2str(i),'.txt'),'w');
fprintf(fileID,'%6s %12s\n','%W[nm]','S[counts]');
for nn=1:size(A,1)
    fprintf(fileID,'%12.8f %12.8f\n',A(nn,:));
end
fclose(fileID);
clear A
%% export 633 excited data
for i = 1:size(spectra1,2)
    DATA = spectra1(:,i)-background1(:);
    DD(i,:) = DATA;
end
%% to save 633
i = 26
spectra_export = DD(i,:);

figure
plot(lambda,spectra_export)
%
A = [lambda,spectra_export'];
save_folder = 'D:\Martin\Publications\anti-stokes\Figures\Supplementary\09_Error_vs_Wavelength\test\';
fileID = fopen(strcat(save_folder,'633nm_spectra_NR',num2str(i),'.txt'),'w');
fprintf(fileID,'%6s %12s\n','%W[nm]','S[counts]');
for nn=1:size(A,1)
    fprintf(fileID,'%12.8f %12.8f\n',A(nn,:));
end
fclose(fileID);