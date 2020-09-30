function [day_data, day_time, night_data, night_time]=F_separate_data(data,timeseries,days,lat,lon)

% INPUTS: 
% data- column vector
% timeseries - column vector

day_data=[]; night_data=[];
day_time=[]; night_time=[];

daymin=24*60;

for i=1:length(days)
    sun_times = F_suncycle(lat,lon,days(i),2880);
    startt=(i-1)*daymin+1; endt=i*daymin;
    time=timeseries(startt:endt);
    datas=data(startt:endt);

    
    if sun_times(2) > sun_times(1)
    dayi1=find(time > (days(i)+sun_times(1)/24));
    dayi=find(time(dayi1) <  (days(i)+sun_times(2)/24) );
    day_data=[day_data; datas(dayi)]; day_time=[day_time; time(dayi)];
    
    nighti1=find(time < (days(i)+sun_times(1)/24));
    nighti2= find(time > (days(i)+sun_times(2)/24));
    night_data=[night_data; datas(nighti1); datas(nighti2)]; 
    night_time=[night_time; time(nighti1); time(nighti2)];
    
    elseif sun_times(2) < sun_times(1)
    nighti1=find(time > (days(i)+sun_times(2)/24));
    nighti=find(time(nighti1) <  (days(i)+sun_times(1)/24) );
    night_data=[night_data; datas(nighti)]; 
    night_time=[night_time; time(nighti)];
    
    dayi1=find(time < (days(i)+sun_times(2)/24));
    dayi2= find(time > (days(i)+sun_times(1)/24));
    day_data=[day_data; datas(dayi1); datas(dayi2)]; 
    day_time=[day_time; time(dayi1); time(dayi2)];
    end
end