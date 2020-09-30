function[x_data, y_data, z_data , X_ID, Y_ID, Z_ID, obj, ncid] = read_geomag_netcdf(ncfname, start_index, n_length, flag)
% The code reads the Geomagnetic data in the netCDF file and returns the
% variable. flag = 0 read only, flag = 1 read/write

S = dir(ncfname);

if ~isempty(S),
    if flag == 0
        ncid = netcdf.open(ncfname,'NOWRITE');
    else
        ncid = netcdf.open(ncfname,'NC_WRITE');
        
        
    end;
    
    if n_length == 0,
        start_index = 0;
        dimid = netcdf.inqDimID(ncid,'data_time_axis_dim')
        [~,n_length] = netcdf.inqDim(ncid,dimid);
        
    end
    
    X_ID = netcdf.inqVarID(ncid,'Magnetic_Field_X');
    Y_ID = netcdf.inqVarID(ncid,'Magnetic_Field_Y');
    Z_ID = netcdf.inqVarID(ncid,'Magnetic_Field_Z');
    
    
    x_data = netcdf.getVar(ncid, X_ID, start_index, n_length);
    x_data = double(x_data)/10;
    x_data(x_data==99999.9) = NaN;
    
    y_data = netcdf.getVar(ncid, Y_ID, start_index, n_length);
    y_data = double(y_data)/10;
    y_data(y_data==99999.9) = NaN;
    
    z_data = netcdf.getVar(ncid, Z_ID, start_index, n_length);
    z_data = double(z_data)/10;
    z_data(z_data==99999.9) = NaN;
    
    [obj] = read_netCDF_header(ncid,1);
    if flag == 0
        
        netcdf.close(ncid);
    end;
    
else
    
    x_data = NaN;
    z_data = NaN;
    y_data = NaN;
    X_ID = NaN;
    Y_ID = NaN;
    Z_ID = NaN;
    ncid = NaN;
    obj = NaN;
    fprintf('Failed to read the file \n');
end;


%******************************************************************
    function [obj] = read_netCDF_header(ncid,full)
        if nargin < 2
            full = 0;
        end
        
        %   get information about global attributes
        obj.id = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'id');
        obj.station_name = netcdf.getAtt(ncid,...
            netcdf.getConstant('NC_GLOBAL'),'station_name');
        obj.geospatial_lat = netcdf.getAtt(ncid,...
            netcdf.getConstant('NC_GLOBAL'),'geospatial_lat');
        obj.geospatial_lon = netcdf.getAtt(ncid,...
            netcdf.getConstant('NC_GLOBAL'),'geospatial_lon');
        obj.geospatial_alt = netcdf.getAtt(ncid,...
            netcdf.getConstant('NC_GLOBAL'),'geospatial_alt');
        obj.time_coverage_start = netcdf.getAtt(ncid,...
            netcdf.getConstant('NC_GLOBAL'),'time_coverage_start_unix_sec_since_1970_01_01');
        obj.time_coverage_end = netcdf.getAtt(ncid,...
            netcdf.getConstant('NC_GLOBAL'),'time_coverage_end_unix_sec_since_1970_01_01');
        obj.time_coverage_duration = netcdf.getAtt(ncid,...
            netcdf.getConstant('NC_GLOBAL'),'time_coverage_resolution_sec');
        
        %             obj.time_coverage_start = netcdf.getAtt(ncid,...
        %                 netcdf.getConstant('NC_GLOBAL'),'time_coverage_start');
        %             obj.time_coverage_end = netcdf.getAtt(ncid,...
        %                 netcdf.getConstant('NC_GLOBAL'),'time_coverage_end');
        %             obj.time_coverage_duration = netcdf.getAtt(ncid,...
        %                 netcdf.getConstant('NC_GLOBAL'),'time_coverage_resolution');
        
        
        if full
            obj.history = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'history');
            obj.date_created = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'date_created');
            obj.creator_name = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'creator_name');
            obj.creator_url = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'creator_url');
            obj.institution = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'institution');
            obj.project = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'project');
            obj.processing_level = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'processing_level');
            obj.station_country = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'station_country');
            obj.station_institution = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'station_institution');
            obj.station_institution_url = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'station_institution_url');
            obj.station_institution_email = netcdf.getAtt(ncid,...
                netcdf.getConstant('NC_GLOBAL'),'station_institution_email');
            
            
            
            
        end
    end
end

