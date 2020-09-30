function [All]=obs_RLSI_timewindows(data,t,lat,lon, startd, endd)
%  function [All, Day, Night]=obs_RLSI_timewindows(data,t,days,lat,lon, startd, endd)
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

time=[]; t1=[];
d=[]; d1=[];
 day=[];

for i=1:n
    day_ind1=find(days==floor(startd(i)));
    day_ind2=find(days==floor(endd(i)));
    day=[day; days(day_ind1:day_ind2)];

    ind1=find(t==startd(i));
    ind2=find(t==endd(i));
    
    temp=(data(ind1:ind2));
    inv=isnan(temp);
    temp1=temp(~inv);
    d=[d; detrend(temp1)];
     d1=[d1; detrend(temp)]; % for separating days/nights
    timetemp=t(ind1:ind2);
    timetemp1=timetemp(~inv);
    time=[time; timetemp1];
    t1=[t1; timetemp]; % for separating days/nights
end
 % periods=Period';
period=[4 4.8 6 8 11.9672 12 12.4206 12.6583 23.9345 24]; %  24.0659 25.8193 26.8683
% period = [ 6 8 12  12.421 23.934472 24];  

 L = isnan(d); % assign which component to process
 
sdata = d(~L); % In this case (NGK), data gaps were in the first and last part.
t = time(~L); % so just delete is OK. Interpolate if data gaps are in between

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
 
%     % Separate day and night
%     [day_data, day_time, night_data, night_time]=F_separate_hr_data(d1,t1,day,lat,lon);
%     
%     % remove NaNs
%     testd=isnan(day_data);
%     day_data=day_data(~testd);
%     day_time=day_time(~testd);
%     testn=isnan(night_data);
%     night_data=night_data(~testn);
%     night_time=night_time(~testn);
%     
%     % For day only
%     disp('start day only')
%     x=[cos(2*pi*day_time/(period(1)/24)) sin(2*pi*day_time/(period(1)/24))];
% 
%     for j=2:length(period)
%         x=[x cos(2*pi*day_time/(period(j)/24)) sin(2*pi*day_time/(period(j)/24))]; 
%     end
% 
%     [spectra_rob,stats] = robustfit(x, day_data);
%     clear x
%     amplitudes = abs(complex(spectra_rob(2:2:end),spectra_rob(3:2:end))); 
%     Day=struct('spectra_rob',spectra_rob, 'stats',stats, 'amplitudes', amplitudes);
%  
%     % For night only
%     disp('start night only')
%     x=[cos(2*pi*night_time/(period(1)/24)) sin(2*pi*night_time/(period(1)/24))];
% 
%     for j=2:length(period)
%         x=[x cos(2*pi*night_time/(period(j)/24)) sin(2*pi*night_time/(period(j)/24))]; 
%     end
% 
%     [spectra_rob,stats] = robustfit(x, night_data);
%     clear x
%     amplitudes = abs(complex(spectra_rob(2:2:end),spectra_rob(3:2:end))); 
%     Night=struct('spectra_rob',spectra_rob, 'stats',stats, 'amplitudes', amplitudes);
 