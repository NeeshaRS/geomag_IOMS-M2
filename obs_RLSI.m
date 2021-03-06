
 function [All]=obs_RLSI(data,t,lat,lon, startd, endd)
% Inputs:
% data - the data series to be inverted
% t - the corresponding time points
% days - the corresponding days
% lat- the station's latitude
% lon- the station's longitude
% startd - start date
% endd - end date

% Outputs:
% All: amplitudes and statistical data for the RLSI on all the data
% Day: amplitudes and statistical data for the RLSI on day-time only data
% Night: amplitudes and statistical data for the RLSI on night-time only data
% ratio: the station's ratio of NaNs to data for the time period

% load Tidal_Periods_hr.mat;

n=length(startd);

time=[];
d=[];

for i=1:n

    ind1=find(t==startd(i));
    ind2=find(t==endd(i));
    
    temp=(data(ind1:ind2));
    inv=isnan(temp);
    temp1=temp(~inv);
    d=[d; detrend(temp1)];
    timetemp=t(ind1:ind2);
    timetemp1=timetemp(~inv);
    time=[time; timetemp1];
end
 % periods=Period';
period=[4 4.8 6 8 11.9672 12 12.4206 12.6583 23.9345 24]; %  24.0659 25.8193 26.8683
%  period = [4 6 8 12  12.4206 23.934472 24];  

 L = isnan(d); % assign which component to process
 
sdata = d(~L);
t = time(~L); 

    % For day and night
    % Model
    m=[cos(2*pi*t(:)/(period(1)/24) ) sin(2*pi*t(:)/(period(1)/24) )];
 
    for j=2:length(period)
    m=[m cos(2*pi*t(:)/(period(j)/24) ) sin(2*pi*t(:)/(period(j)/24) )]; 
    end
 
    % Now invert
    [spectra_rob, stats] = robustfit(m,detrend(sdata));

    clear x
    amplitudes = abs(complex(spectra_rob(2:2:end),spectra_rob(3:2:end))); 
    All=struct('spectra_rob',spectra_rob, 'stats',stats, 'amplitudes', amplitudes);

 