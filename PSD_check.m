%% Make power spectra
% load and read the observatory mat files
close all; clc; clear all;
load Obs_Time.mat 
% localdirectory = '/home/nschnepf/Ground_Observatories/';
addpath /home/nschnepf/Ground_Observatories/raw_data/
S = dir('/home/nschnepf/Ground_Observatories/raw_data/*raw_data.mat');

for i=1:length(S);
    load(S(i).name)
    name_index(i,1:3)=S(i).name(1:3);
end

whos
%%
names=['API'; 'ASP'; 'GNA'; 'GUA'; 'KAK'; 'KNY'; 'MMB'; 'NEW'; 'PPT'; 'TAM'; 'WNG'];
index=[8 12 57 59 77 82 106 109 121 135 153];

clc
tic
for i=index;
    load(S(i).name)
    
    % startd= datenum('January 1 1995'); % start time for the RLSI
    % endd= datenum('December 31 1996 23:59:00'); % end time for the RLSI
    
    % Solar minimum 1
     startd= datenum('January 1 1995'); % start time for the RLSI
     endd= datenum('December 31 1998 23:59:00'); % end time for the RLSI
    
    % Solar maximum
    % startd= datenum('January 1 2000'); % start time for the RLSI
    % endd= datenum('December 31 2003 23:59:00'); % end time for the RLSI
    
    % Solar minimum 2
    % startd= datenum('January 1 2007'); % start time for the RLSI
    % endd= datenum('December 31 2010 23:59:00'); % end time for the RLSI
    
     tic
     name=S(i).name(1:3);
     F_obs_psd(F,t, startd, endd,name);
     toc
    

end
toc