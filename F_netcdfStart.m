function [] = F_netcdfStart(S1, directory, outnet)

% Inputs:
% S1: specifies the directory that has compressed netcdf files. 
% Ex. S1 = dir('/data/obs_mag_data/update_2012/absolute/*.gz');
% directory: a string of the directory of the compressed netcdf files.
% Ex. directory = '/data/obs_mag_data/update_2012/absolute/';
% outnet: a string of the directory for the output unzipped netcdf files
% Ex. outnet='/home/nschnepf/Ground_Observatories/';

% Unzip all the files
for i = 1:length(S1)
    gunzip([directory S1(i).name],outnet);
end

% load and read the unzipped observatory netcdf files
S = dir([outnet '*.nc']);
coords=zeros(length(S),2);

for i=1:length(S)
    ncid=netcdf.open(S(i).name,'NOWRITE');
    [X, Y, Z , X_ID, Y_ID, Z_ID, obj, ncid] = read_geomag_netcdf(S(i).name, 0, n_length, 0);
    F=(X.^2 + Y.^2 + Z.^2).^.5;
    % deteremine the station coordinates lat, lon
    lat=double(obj.geospatial_lat); lon=double(obj.geospatial_lon);
    if length(lat)>1
        lat=str2double(obj.geospatial_lat); 
        lon=str2double(obj.geospatial_lon);
    end
    fname=strcat(S(i).name(1:3),'_raw_data.mat')
    save(fname,'X','Y', 'Z', 'F', 'lat', 'lon')
    
    coords(i,1)=lat;
    coords(i,2)=lon;
end

save('stations_coords.mat','coords')
disp('stations_coods.mat saved')