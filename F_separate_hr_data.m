function [day_data, day_time, night_data, night_time]=F_separate_hr_data(data,timeseries,days,lat,lon)

day_data=[]; night_data=[];
day_time=[]; night_time=[];

dayhr=24; % for hourly data

for i=1:length(days)
    sun_times = F_suncycle(lat,lon,days(i),2880);
    startt=(i-1)*dayhr+1; endt=i*dayhr;
    time=timeseries(startt:endt);
    datas=data(startt:endt);

    
    if sun_times(2) > sun_times(1)
    dayi1=find(time > (days(i)+sun_times(1)/24) & time  <  (days(i)+sun_times(2)/24));
    day_data=[day_data; datas(dayi1)]; day_time=[day_time; time(dayi1)];
    
    nighti1=find(time < (days(i)+sun_times(1)/24) | time > (days(i)+sun_times(2)/24));
    night_data=[night_data; datas(nighti1)]; 
    night_time=[night_time; time(nighti1)];
    
    elseif sun_times(2) < sun_times(1)
    nighti=find(time > (days(i)+sun_times(2)/24) & time <  (days(i)+sun_times(1)/24));
    night_data=[night_data; datas(nighti)]; 
    night_time=[night_time; time(nighti)];
    
    dayi1=find(time < (days(i)+sun_times(2)/24) | time > (days(i)+sun_times(1)/24));
    day_data=[day_data; datas(dayi1)]; 
    day_time=[day_time; time(dayi1)];
    end
end