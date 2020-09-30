
function [p] = F_compare_PSD_RLSI(startd, endd, wl, data, t, name, component)
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
    
    check=isempty(temp1);
    if check == 0
    d=[d; detrend(temp1)];
    timetemp=t(ind1:ind2);
    timetemp=timetemp(~inv);
    time=[time; timetemp];
    
    temp(inv)=0; % change NaNs to zeros
    fd=temp+fd;
    
    end
end
fd=fd/n;

L = isnan(d); 
sdata = d(~L); 
t = time(~L); 

check=isempty(L);
    if check == 1
        disp('not enough data points')
        p=[1 1];
       return
    end

% Recover amplitude and phase using RLSI fitting
periods = [ 6 8 12  12.421 23.934472 24].*3600; % Try adding and removing tides
period =[4 4.8 6 8 11.9672 12 12.4206 12.6583 23.9345 24]; 

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
hold on; % tidal modes for the RLSI
line([4 4],[2 14],'LineWidth',.5,'Color','r');
line([4.8 4.8],[2 14],'LineWidth',.5,'Color','r');
line([6 6],[2 14],'LineWidth',.5,'Color','r');
line([8 8],[2 14],'LineWidth',.5,'Color','r');
line([11.967 11.967],[2 14],'LineWidth',.5,'Color','r');
line([12 12],[2 14],'LineWidth',.5,'Color','r');
line([12.421 12.421],[2 14],'LineWidth',.5,'Color','r');
line([12.658 12.658],[2 14],'LineWidth',.5,'Color','r');
line([23.934 23.934],[2 14],'LineWidth',.5,'Color','r');
line([24 24],[2 14],'LineWidth',.5,'Color','r');

plot(1./(fr.*3600),abs(half_spectra),'LineWidth',2,'Color','k');
set(gca,'FontSize',15)
xlabel('Period (hr)'); ylabel('FFT Amplitude (nT)')
xlim([2 25])

year1=datestr(startd(1)); year1=year1(8:11);
year2=datestr(endd(end)); year2=year2(8:11);
figname=[name '_' year1 '-' year2 '_' num2str(wl) '_' component '.png'];
print(figname,'-dpng')

% Get FFT spectral lines closest to the tidal periods 
for i = 1:length(periods),
    fft_index(i) = findnearest(periods(i),1./fr,0);
end;

% Now print the results to a file
fname=[name '_' year1 '-' year2 '_' num2str(wl) '_' component '.txt'];
fid=fopen(fname, 'wt');

for i = 1:length(periods),
    %fprintf('%5.4f  %5.4f %5.4f %5.4f\n', periods(i)/3600, abs(complex(a(2*i),a(2*i+1))), ...
    %     stats.p(2*i),stats.p(2*i+1));
    fprintf(fid,'%5.4f %2.4f%5.2f %5.4f %5.4f %5.4f\n', periods(i)/3600, 1/(3600*fr(fft_index(i))),abs(half_spectra(fft_index(i))), abs(complex(a(2*i),a(2*i+1))), ...
        stats.p(2*i),stats.p(2*i+1));
end;

p=[stats.p(2*7),stats.p(2*7+1)];
