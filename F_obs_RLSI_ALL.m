function [All]=F_obs_RLSI_ALL(data,t)
load Tidal_Periods_hr.mat;

% periods=[4 4.8 6 8 11.967236 12 12.42060122 12.6583 ... 
%     23.934472 24 24.0659 25.8913 25.8956 4382.91 8765.81 nm];
% period_string=['S6';'S5';'S4';'S3';'K2';'S2';'M2';'N2';...
%     'K1';'S1';'P1';'O1';'Q1';'HY';'Yr';'NM'];

periods=Period';

% For day and night
disp('start day and night')
x=[cos(2*pi*t/(periods(1)/24)) sin(2*pi*t/(periods(1)/24))];

for j=2:length(periods)
   x=[x cos(2*pi*t/(periods(j)/24)) sin(2*pi*t/(periods(j)/24))]; 
end

 [spectra_rob,stats] = robustfit(x, data);
 clear x
 amplitudes = abs(complex(spectra_rob(2:2:end),spectra_rob(3:2:end))); 
 All=struct('spectra_rob',spectra_rob, 'stats',stats, 'amplitudes', amplitudes);