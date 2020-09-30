function [day_data, day_time]=F_separate_data(data,timeseries,days,lat,lon)

for i=1:length(days)
    sun_times = F_suncycle(lat,lon,days(i),2880);
end