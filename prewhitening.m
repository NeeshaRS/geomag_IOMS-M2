%% load in data from IRT
clc; clear all; close all;

% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/ASP_1995_2015.nc'
% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/BOU_1995_2015.nc'
 ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/BDV_1995_2015.nc'
% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/NGK_1995_2015.nc'
% ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/NVS_1995_2015.nc'
%  ncfname = '/data/nschnepf/Geomag_Observatories/Hourly_Data/TSU_1995_2015.nc'

t=datenum(1995,1,1,0,0,30):(1/24):datenum(2015,12,31,23,59,30); t=t';
n1 = length(t);% length of time series in the netCDF file
[x_data_h, y_data_h, z_data_h , X_ID, Y_ID, Z_ID, obj] = read_geomag_netcdf(ncfname, 0, n1, 0);% Get data

X=x_data_h;
Y=y_data_h;
Z=z_data_h;

clear *data*

disp('data loaded')
%% 1: Fourier Transform
L=isnan(Z);
data=Z(~L);
time=t(~L);

fs=1/3600; % 1 hour sampling
numsamples=length(data);

data_fft=fft(detrend(data)); % Fourier Transform
fr=fftfrq(numsamples,fs); % Get corresponding frequencies

% data_fft=data_fft(2:numsamples/2+1)/(numsamples/2); % scaling
plot(1./(fr.*3600),data_fft(2:numsamples/2+1)/(numsamples/2)); xlim([0 25])

phase=angle(data_fft);

disp('1: fourier transform complete')
%% 2: SMOOTH POWER SPECTRA
NW=5.5; %2, 5/2, 3, 7/2, 4 are the typical values
NFFT=numsamples;
[Pxx, w] = pmtm(data_fft,NW,NFFT);
Pxx=Pxx';
w=w';
wHz=w./(2*pi*3600); % convert to Hz

Pxx_mag=sqrt(ifftshift(Pxx)); % the smoothed magnitude
% Pxx_mag=sqrt(Pxx); % the smoothed magnitude

whos
disp('2: smooth spectra made')
%% 3: Whiten
temp=Pxx_mag/max(Pxx_mag); % Scale the smooth spectra
white_mag=abs(data_fft)'./temp;
white_freqDom=white_mag'.*exp(i*phase); % pre-whitened signal in frequency domain
% white_freqDom=white_freqDom/(numsamples/2); % scale it
white_timeDom=ifft(white_freqDom);
whos

plot(1./(fr.*3600),abs(white_freqDom(2:numsamples/2+1))); xlim([5 15]);
half_spectra=white_freqDom(2:numsamples/2+1);

disp('3: whitened')
%% 4: compare RLSI M2 amplitude for raw & whitened data
periods = [ 6 8 12  12.421 23.934472 24].*3600; %
t=(1:numsamples)/fs;

% Use RLSI on raw data
m=[cos(2*pi*t(:)/(periods(1)) ) sin(2*pi*t(:)/(periods(1)) )];
for j=2:length(periods)
   m=[m cos(2*pi*t(:)/(periods(j)) ) sin(2*pi*t(:)/(periods(j)) )]; 
end
[a, stats] = robustfit(m,detrend(data));

% Use RLSI on pre-whitened data
t=t(20:end-19);
m1=[cos(2*pi*t(:)/(periods(1)) ) sin(2*pi*t(:)/(periods(1)) )];
for j=2:length(periods)
   m1=[m1 cos(2*pi*t(:)/(periods(j)) ) sin(2*pi*t(:)/(periods(j)) )]; 
end
[a1, stats1] = robustfit(m1,(abs(white_timeDom(20:end-19))));
 
% Now print the results
for i = 1:length(periods),
    fprintf('%5.4f %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n', periods(i)/3600, abs(complex(a(2*i),a(2*i+1))),...
        stats.p(2*i), stats.p(2*i+1), abs(complex(a1(2*i),a1(2*i+1))),...
        stats1.p(2*i), stats1.p(2*i+1));
end;