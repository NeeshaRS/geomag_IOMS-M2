
function [] = compare_PSD_RLSI(startd, endd, wl, data, t, name, component)
% Inputs:
% startd: the start dates for each time window
% endd: the end dates for each time window
% wl: the length of each window
% data: the data series
% t: the corresponding time series
% name: the name of the station [string]
% component: the magnetic field component being looked at [string]

n=length(startd);
time=[];
d=[];

fd=zeros(wl,1);

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

L = isnan(d); 
sdata = d(~L); 
t = time(~L); 

% Recover amplitude and phase using RLSI fitting
periods = [ 6 8 12  12.421 23.934472 24].*3600; % Try adding and removing tides
period = [ 6 8 12  12.421 23.934472 24]; 

m=[cos(2*pi*t(:)/(period(1)/24) ) sin(2*pi*t(:)/(period(1)/24) )];
for j=2:length(period)
   m=[m cos(2*pi*t(:)/(period(j)/24) ) sin(2*pi*t(:)/(period(j)/24) )]; 
end

[a, stats] = robustfit(m,detrend(sdata));

% Recover amplitude using fft
fs = 1/3600; % sampling every 1 hour
numsamples = length(fd);
tf = (1:numsamples)/fs; %(relative time Axis in seconds)

c = fft(fd); % Fourier transform
fr = fftfrq(numsamples, fs); % Get frequency for FFT 
half_spectra = c(2:numsamples/2+1)/(length(fd)/2); % get first half of conjugate spectra and scale it with N samples

h=figure(1);
plot(1./(fr.*3600),abs(half_spectra),'LineWidth',2);
set(gca,'FontSize',15)
xlabel('Period (hr)'); ylabel('FFT Amplitude')
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
