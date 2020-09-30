%% The data spans from the start of 1995 to the end of 2012
clc; clear all; close all; 

start_t=datenum('January 1 1995');
end_t=datenum('January 1 2013');
ts=1/(24*60);
t=start_t:ts:end_t-ts;
day=start_t:end_t-ts;
save Obs_Time.mat t day
n_length=length(t);

%% Get tidal amplitudes
% load and read the observatory mat files
close all; clc; clear all;
load Obs_Time.mat 
% localdirectory = '/home/nschnepf/Ground_Observatories/';
addpath /home/nschnepf/Ground_Observatories/raw_data/
S = dir('/home/nschnepf/Ground_Observatories/raw_data/*raw_data.mat');

% get the tidal amplitudes
for i=1:length(S);
    tic
    load(S(i).name)
    
    % r= 2/3;
     r= .5;
    
    % Solar minimum 1
    % startd= datenum('January 1 1997'); % start time for the RLSI
    % endd= datenum('December 31 1998 23:59:00'); % end time for the RLSI
    % endd= datenum('December 31 1998 23:59:00'); % end time for the RLSI
    
    % Solar maximum
    % startd= datenum('January 1 2000'); % start time for the RLSI
    % endd= datenum('December 31 2003 23:59:00'); % end time for the RLSI
    
    % Solar minimum 2
    % startd= datenum('January 1 2007'); % start time for the RLSI
    % endd= datenum('December 31 2010 23:59:00'); % end time for the RLSI
    
    [All, Day, Night, dratio]=obs_RLSI(F,t',day,lat,lon, startd,endd,r);
    fname=strcat(S(i).name(1:3),'_tidal_analysis_1997-1998.mat')
    save(fname,'All','Day','Night','dratio')
    toc
    
    clear all;
    load Obs_Time.mat 
    S = dir('/home/nschnepf/Ground_Observatories/raw_data/*raw_data.mat');
end

%% Get tidal amplitudes for summer days only
% load and read the observatory mat files
close all; clc; clear all;
addpath /data/nschnepf/Geomag_Observatories/Hourly_Data/
S = dir('/data/nschnepf/Geomag_Observatories/Hourly_Data/*.nc'); % hourly data

load good_stations_SM_2009_Y.mat
n=length(good_stations.index);
index=zeros(n,1);
for i=1:n
   index(i)=good_stations.index{i};
end
save index_2009.mat index

tic
% get the tidal amplitudes
    
for i=1:n % process the good stations
    tic

    t=datenum(1995,1,1,0,0,30): (1/(24)): datenum(2015,12,31,23,59,30); t=t';
    day=datenum(1995,1,1,0,0,0): (1): datenum(2015,12,31,23,59,0); day=day';
    n1 = length(t);% length of time series in the netCDF file
    path='/data/nschnepf/Geomag_Observatories/Hourly_Data/';
    
    load index_2009.mat
%     v=i; % for individual tests
    v=index(i);
    name=[path S(v).name]
    [X, Y, Z , X_ID, Y_ID, Z_ID, obj] = read_geomag_netcdf(name, 0, n1, 0);
    F=(Y.^2+X.^2+Z.^2).^.5;
    
    lat=obj.geospatial_lat
    lon=obj.geospatial_lon
    
startd=[datenum('May 30 2009 00:00:30')];
endd=startd+90;
    
    if abs(lat) < 56 % ignore polar stations
        [All]=F_obs_RLSI_timewindows(F,t,lat,lon, startd,endd);
            fname=strcat(S(v).name(1:3),'_F_10m_TA_07-09.mat')
            save(fname,'All','lat','lon') %
        toc
    end
    clear all; 
    S = dir('/data/nschnepf/Geomag_Observatories/Hourly_Data/*.nc');
    % S = dir('/home/nschnepf/Ground_Observatories/raw_data/*raw_data.mat');
end
toc
disp('done')
%% Compare the PSD and RLSI to compare spectra & p-values to determine if
% the station is worth pursuing

close all; clc; clear all;
addpath /data/nschnepf/Geomag_Observatories/Hourly_Data/
S = dir('/data/nschnepf/Geomag_Observatories/Hourly_Data/*.nc'); % hourly data

startd=[datenum('May 30 2009 00:00:30')];
endd=startd+90;

wl=length(startd:1/24:endd); % number of hourly data points during each time window

t=datenum(1995,1,1,0,0,30): (1/(24)): datenum(2015,12,31,23,59,30); t=t';
n1 = length(t);% length of time series in the netCDF file

tic

bad_p_stations=[];
k=1; % counter
for i=21; %1:length(S)
    
    path='/data/nschnepf/Geomag_Observatories/Hourly_Data/';
    name=[path S(i).name];
    [x_data_h, y_data_h, z_data_h , X_ID, Y_ID, Z_ID, obj] = read_geomag_netcdf(name, 0, n1, 0);
    F=(y_data_h.^2+x_data_h.^2+z_data_h.^2).^.5;
    
    lat=obj.geospatial_lat
    lon=obj.geospatial_lon
    
    if abs(lat) < 56 % ignore polar stations
        name=S(i).name(1:3)
        component='Y'
        p= F_compare_PSD_RLSI(startd, endd, wl, y_data_h, t, name,component);
        if (p(1)^2 + p(2)^2)^.5 >= 0.05
            bad_p_stations{k}=S(i).name(1:3);
            disp(bad_p_stations{k})
            k=k+1;
        end
    end
end
%% Add stations with bad spectra
names=['API'; 'ASC'; 'ASP'; 'FRD'; 'FRN'; 'GCK'; 'HBK'; 'IZN'; 'KOU'; 'MBO'; ...
    'MEA'; 'NMP'; 'OTT'; 'PAF'; 'PAG'; 'PND'; 'PPT'; 'SFS'; 'SHU'; 'SJG'; 'STJ'; ...
    'THY'; 'TRW'; 'VIC'; 'VSS'];
for i=1:length(names)
    bad_p_stations{k}=names(i,:);
    k=k+1;
end
toc
save stations_with_bad_p_vals_SM_2009_Y.mat bad_p_stations
disp('done')

%% Determine stations to pursue

close all; clc; clear all;
addpath /data/nschnepf/Geomag_Observatories/Hourly_Data/
S = dir('/data/nschnepf/Geomag_Observatories/Hourly_Data/*.nc'); % hourly data

load stations_with_bad_p_vals_SM_2009_Y.mat

good_stations=[];
k=1; % counter

for i=1:length(S)
    name=S(i).name(1:3)
    all_stations{i}=name;
    test=0;
   for j=1:length(bad_p_stations)
       test= test+strcmp(name,bad_p_stations{j});
   end
   
   if test==0
       good_stations.station{k}=name;
       good_stations.index{k}=i;
       k=k+1;
   end
   
end
save good_stations_SM_2009_Y.mat good_stations all_stations
disp('done')

%% Make power spectra
% load and read the observatory mat files
close all; clc; clear all;
% load Obs_Time.mat 
% localdirectory = '/home/nschnepf/Ground_Observatories/';
% addpath /home/nschnepf/Ground_Observatories/raw_data/
% S = dir('/home/nschnepf/Ground_Observatories/raw_data/*raw_data.mat');

addpath /data/nschnepf/Geomag_Observatories/Hourly_Data/
S = dir('/data/nschnepf/Geomag_Observatories/Hourly_Data/*.nc'); % hourly data

%% get the tidal amplitudes
tic
for i=21; %1:length(S);
    load(S(i).name)
    
    startd=[datenum('May 30 2009 00:00:30')];
    endd=startd+90;
    
    % Solar minimum 1
    % startd= datenum('January 1 1995'); % start time for the RLSI
    % endd= datenum('December 31 1998 23:59:00'); % end time for the RLSI
    
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