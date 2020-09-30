%% Load in observatories
% load and read the observatory mat files
close all; clc; clear all;

% addpath 2009_Summer/
% S= dir('/home/nschnepf/Ground_Observatories/2009_Summer/*F_8m_TA_2009.mat');

S= dir('/home/nschnepf/Ground_Observatories/2009_Summer/*_Y_10m_TA_2009.mat');
%% stats
numstats=1+1+1+1+1; % amplitude + standard error + p values + lat + lon
obs_ALL_stats=zeros(length(S),numstats);
% obs_DAY_stats=zeros(length(S),numstats);
% obs_NIGHT_stats=zeros(length(S),numstats);

tic
final_stations=[];
for k=1:length(S)
    load(S(k).name)
        
        obs_ALL_stats(k,1)=All.amplitudes(7);
        obs_ALL_stats(k,2)=(All.stats.se(15)^2+All.stats.se(16)^2)^.5;
        obs_ALL_stats(k,3)=(All.stats.p(15)^2+All.stats.p(16)^2)^.5;
        obs_ALL_stats(k,4)=lat;
        obs_ALL_stats(k,5)=lon;
        
%         obs_DAY_stats(k,1)=Day.amplitudes(4);
%         obs_DAY_stats(k,2)=(Day.stats.se(8)^2+Day.stats.se(9)^2)^.5;
%         obs_DAY_stats(k,3)=(Day.stats.p(8)^2+Day.stats.p(9)^2)^.5;
%         obs_DAY_stats(k,4)=lat;
%         obs_DAY_stats(k,5)=lon;
        
%         obs_NIGHT_stats(k,1)=Night.amplitudes(4);
%         obs_NIGHT_stats(k,2)=(Night.stats.se(8)^2+Night.stats.se(9)^2)^.5;
%         obs_NIGHT_stats(k,3)=(Night.stats.p(8)^2+Night.stats.p(9)^2)^.5;
%         obs_NIGHT_stats(k,4)=lat;
%         obs_NIGHT_stats(k,5)=lon;
        
        final_stations{k}=S(k).name(1:3);
end
toc

save M2_Y_10m_obs_stats_p_2009.mat obs_ALL_stats final_stations

disp('saved')