clc; clear all; close all;
% Script to get M2 signal amplitude from fft and inversion
% Get hourly data
% ncfname = '/Users/manojnair/data/obs_mag_data/BGS4/NGK_1995_2015.nc'

%  ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/ASP_1995_2015.nc'
% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/BOU_1995_2015.nc'
ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/HON_1995_2015.nc'
% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/IRT_1995_2015.nc'
% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/NGK_1995_2015.nc'
% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/NVS_1995_2015.nc'
%  ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/TSU_1995_2015.nc'

t=datenum(1995,1,1,0,0,30): (1/(24)): datenum(2015,12,31,23,59,30); t=t';
n1 = length(t);% length of time series in the netCDF file
fday_h = datenum(1995,1,1,0,0,30): (1/(24)): datenum(2015,12,31,23,59,30);% Get absolute time axis
[x_data_h, y_data_h, z_data_h , X_ID, Y_ID, Z_ID, obj] = read_geomag_netcdf(ncfname, 0, n1, 0);% Get data
plot(fday_h, x_data_h,'r'); % Plot data to check any major issue with data (like baseline0
 
whos
%% Treat data gaps
 
L = isnan(y_data_h); % assign which component to process
 
sdata = y_data_h(~L); % In this case (NGK), data gaps were in the first and last part.
fday_n = fday_h(~L); % so just delete is OK. Interpolate if data gaps are in between
 
 
% Make a smaller set of  time series  
% 
% sdata_1 = sdata(1:length(sdata)/20);
% fday_n_1 = fday_n(1:length(sdata)/20);
% 
% sdata = sdata_1;
% fday_n = fday_n_1;
 
% If required, plot Power spectra
% fs = 1/3600; % sampling every 1 hour
% nfft = 1024*10;
% [Pyy,F] = pwelch(sdata, hamming(nfft),1,nfft,fs);
% loglog(1./(F*3600),abs(Pyy),'k', 'LineWidth', 2);
 
 
% Invert
fs = 1/3600; % sampling every 1 hour
numsamples = length(sdata);
t = (1:numsamples)/fs; %(relative time Axis in seconds)
 
%2) Recover amplitude and phase using fitting
 
periods = [ 6 8 12  12.421 23.934472 24].*3600; % Try adding and removing tides
 
% Model
m=[cos(2*pi*t(:)/(periods(1)) ) sin(2*pi*t(:)/(periods(1)) )];
 
for j=2:length(periods)
   m=[m cos(2*pi*t(:)/(periods(j)) ) sin(2*pi*t(:)/(periods(j)) )]; 
end
 
% Now invert
[a, stats] = robustfit(m,detrend(sdata));
 
 
% Recover amplitude using fft
 
c = fft(detrend(sdata)); % Fourier transform
fr = fftfrq(numsamples, fs); % Get frequency for FFT 
half_spectra = c(2:numsamples/2+1)/(length(sdata)/2); % get first half of conjugate spectra and scale it with N samples
 
plot(1./(fr.*3600),abs(half_spectra)); % Plot
 
path(path,'/Users/manojnair/m/extern/'); % Just adding path where findnearest.m resides
 
% Get FFT spectral lines closest to the tidal periods 
for i = 1:length(periods),
    fft_index(i) = findnearest(periods(i),1./fr,0);
end;
 
% Now print the results
for i = 1:length(periods),
    fprintf('%5.4f %2.4f%5.2f %5.4f %5.4f %5.4f\n', periods(i)/3600, 1/(3600*fr(fft_index(i))),abs(half_spectra(fft_index(i))), abs(complex(a(2*i),a(2*i+1))), ...
        stats.p(2*i),stats.p(2*i+1));
end;
 
%% using gapped data
startd=[datenum('June 30 2007 00:00:30') 
    datenum('June 30 2008 00:00:30') 
    datenum('June 30 2009 00:00:30')
    datenum('June 30 2010 00:00:30')
    datenum('June 30 2011 00:00:30')];

    endd=[datenum('July 28 2007 23:00:30') 
    datenum('July 28 2008 23:00:30') 
    datenum('July 28 2009 23:00:30')
    datenum('July 28 2010 23:00:30')
    datenum('July 28 2011 23:00:30')];

% startd=[datenum('June 1 2007 00:00:30') 
%     datenum('June 1 2008 00:00:30') 
%     datenum('June 1 2009 00:00:30')
%     datenum('June 1 2010 00:00:30')
%     datenum('June 1 2011 00:00:30')];
% 
%     endd=[datenum('August 31 2007 23:00:30') 
%     datenum('August 31 2008 23:00:30') 
%     datenum('August 31 2009 23:00:30')
%     datenum('August 31 2010 23:00:30')
%     datenum('August 31 2011 23:00:30')];

%   startd=[datenum('January 1 2007 00:00:30')];
%   endd=[datenum('December 31 2011 23:00:30')];

n=length(startd);
time=[];
d=[];
% fd=zeros(2208,1); % for summer timespan
 fd=zeros(696,1); % for 29 day time span
% fd=zeros(43824,1); % for the 5 year period

% data= (y_data_h.^2+x_data_h.^2+z_data_h.^2).^.5;
 data= y_data_h; % look at Y component

for i=1:n
    ind1=find(t==startd(i));
    ind2=find(t==endd(i));

    temp=(data(ind1:ind2));
    inv=isnan(temp);
    temp1=temp(~inv);
    d=[d; detrend(temp1)];
    timetemp=t(ind1:ind2);
    timetemp=timetemp(~inv);
    time=[time; timetemp];
    
    temp(inv)=0; % change NaNs to zeros
    fd=temp+fd;
 
end
fd=fd/n;
%
L = isnan(d); % assign which component to process
 
sdata = d(~L); % In this case (NGK), data gaps were in the first and last part.
t = time(~L); % so just delete is OK. Interpolate if data gaps are in between

% Invert
fs = 1/3600; % sampling every 1 hour
numsamples = length(fd);
tf = (1:numsamples)/fs; %(relative time Axis in seconds)

%1) Recover amplitude and phase using fitting
 
periods = [ 6 8 12  12.421 23.934472 24].*3600; % Try adding and removing tides
period = [ 6 8 12  12.421 23.934472 24]; 

% Model
m=[cos(2*pi*t(:)/(period(1)/24) ) sin(2*pi*t(:)/(period(1)/24) )];
 
for j=2:length(periods)
   m=[m cos(2*pi*t(:)/(period(j)/24) ) sin(2*pi*t(:)/(period(j)/24) )]; 
end
 
% Now invert
[a, stats] = robustfit(m,detrend(sdata));
 
% Recover amplitude using fft
 
c = fft(fd); % Fourier transform
fr = fftfrq(numsamples, fs); % Get frequency for FFT 
half_spectra = c(2:numsamples/2+1)/(length(fd)/2); % get first half of conjugate spectra and scale it with N samples

figure(1)
plot(1./(fr.*3600),abs(half_spectra)); % Plot
xlim([10 15])

% Get FFT spectral lines closest to the tidal periods 
for i = 1:length(periods),
    fft_index(i) = findnearest(periods(i),1./fr,0);
end;

% Now print the results
for i = 1:length(periods),
    %fprintf('%5.4f  %5.4f %5.4f %5.4f\n', periods(i)/3600, abs(complex(a(2*i),a(2*i+1))), ...
    %     stats.p(2*i),stats.p(2*i+1));
    fprintf('%5.4f %2.4f%5.2f %5.4f %5.4f %5.4f\n', periods(i)/3600, 1/(3600*fr(fft_index(i))),abs(half_spectra(fft_index(i))), abs(complex(a(2*i),a(2*i+1))), ...
        stats.p(2*i),stats.p(2*i+1));
end;
