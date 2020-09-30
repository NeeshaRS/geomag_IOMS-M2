%% Load in observatories
% load and read the observatory mat files
close all; clc; clear all;

%%
addpath /data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/
S= dir('/data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/*X_10m_TA_2009.mat');
obs_X_phase = zeros(length(S),3);

for k=1:length(S)
    load(S(k).name)
        
        obs_X_phase(k,1)=ata2nd(All.spectra_rob(15), All.spectra_rob(14));
        obs_X_phase(k,2)=lat;
        obs_X_phase(k,3)=lon;
end

save obs_X_phase.mat obs_X_phase
disp('saved')
%%
addpath /data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/
S= dir('/data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/*Y_10m_TA_2009.mat');
obs_Y_phase = zeros(length(S),3);

for k=1:length(S)
    load(S(k).name)
        
        obs_Y_phase(k,1)=atan2d((All.spectra_rob(15)),(All.spectra_rob(14)));
        obs_Y_phase(k,2)=lat;
        obs_Y_phase(k,3)=lon;
end

save obs_Y_phase.mat obs_Y_phase
disp('saved')

plot(obs_Y_phase(:,1),'o')
%%
addpath /data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/
S= dir('/data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/*Z_10m_TA_2009.mat');
obs_Z_phase = zeros(length(S),3);

for k=1:length(S)
    load(S(k).name)
        
        obs_Z_phase(k,1)=atan2d((All.spectra_rob(15)),(All.spectra_rob(14)));
        obs_Z_phase(k,2)=lat;
        obs_Z_phase(k,3)=lon;
end

save obs_Z_phase.mat obs_Z_phase
disp('saved')

plot(obs_Z_phase(:,1),'o')
%%
addpath /data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/
S= dir('/data/nschnepf/Geomag_Observatories/RLSI_Results/2009_Summer/*F_10m_TA_2009.mat');
obs_F_phase = zeros(length(S),3);

for k=1:length(S)
    load(S(k).name)
        
        obs_F_phase(k,1)=atan2d((All.spectra_rob(15)),(All.spectra_rob(14)));
        obs_F_phase(k,2)=lat;
        obs_F_phase(k,3)=lon;
end

save obs_F_phase.mat obs_F_phase
disp('saved')

plot(obs_F_phase(:,1),'o')