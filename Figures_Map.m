%% Initialize
close all; clc; clear all;

% addpath Using_47_modes/1995-1998_Solar_Minimum/
addpath Using_47_modes/2000-2003_Solar_Maximum/
% addpath Using_47_modes/2007-2010_Solar_Minimum/

% S = dir('Using_47_modes/1995-1998_Solar_Minimum/*.mat');
S = dir('Using_47_modes/2000-2003_Solar_Maximum/*.mat');
% S = dir('Using_47_modes/2007-2010_Solar_Minimum/*.mat');
load('M2_Abs_Colormap.mat');
load stations_coords.mat
%% Open background file
open('From_Manoj/M2_Bz_Abs_0_720_1440.fig');
colormap(map1);
caxis([0 10])
%% Plot the M2 day and night amplitudes for each observatory
hold on;

All_p_vals=zeros(length(S),2);
for i = 1:length(S);
    load(S(i).name);
    
    % make sure the station had proper data
    if dratio<2/3
        All_p_vals(i,:)=[All.stats.p(60:61)]';
        m2_amp=All.amplitudes(30);
        
        % Find color to plot
        cidx = floor(63*(m2_amp)/(10)+1);
        if cidx > 64,
            cidx = 64;
        elseif cidx < 1,
            cidx = 1;
        end;
        
        plot(coords(i,2), 90-coords(i,1), '.', 'color',map1(cidx,:),'MarkerSize',50);
    end
    
end
disp('done--> plot of day & night signals')
%% Plot the M2 day amplitudes for each observatory
hold on;

Day_p_vals=zeros(length(S),2);
for i = 1:length(S);
    load(S(i).name);
    
    % make sure the station had proper data
    if dratio<2/3
        Day_p_vals(i,:)=[Day.stats.p(60:61)]';
        m2_amp=Day.amplitudes(30);
        
        % Find color to plot
        cidx = floor(63*(m2_amp)/(10)+1);
        if cidx > 64,
            cidx = 64;
        elseif cidx < 1,
            cidx = 1;
        end;
        
        plot(coords(i,2), 90-coords(i,1), '.', 'color',map1(cidx,:),'MarkerSize',50);
    end
    
end
disp('done--> plot of day signals')

%% Plot the M2 night amplitudes for each observatory
hold on;

Night_p_vals=zeros(length(S),2);
for i = 1:length(S);
    load(S(i).name);
    
    % make sure the station had proper data
    if dratio<2/3
        Night_p_vals(i,:)=[Night.stats.p(60:61)]';
        m2_amp=Night.amplitudes(30);
        
        % Find color to plot
        cidx = floor(63*(m2_amp)/(10)+1);
        if cidx > 64,
            cidx = 64;
        elseif cidx < 1,
            cidx = 1;
        end;
        
        plot(coords(i,2), 90-coords(i,1), '.', 'color',map1(cidx,:),'MarkerSize',50);
    end
    
end
disp('done--> plot of nightsignals')

%% Save p values
save M2_p_values.mat All_p_vals Day_p_vals Night_p_vals
disp('saved')
%% list of stations where data quality was not good enough
clc
gcount=1;
bcount=1;

for i = 1:length(S);
    load(S(i).name);
    
    ratios(i)=dratio;
    % make sure the station had proper data
    if dratio>=2/3
        badlist(bcount,1:3)=S(i).name(1:3);
        bcount=bcount+1;
    elseif dratio<2/3
        goodlist(gcount,1:3)=S(i).name(1:3);
        gratio(gcount)=dratio;
        gcount=gcount+1;
    end
    
end

% save stations_used_list_1995-1998.mat badlist goodlist gratio
% save stations_used_list_2000-2003.mat badlist goodlist gratio
save stations_used_list_2007-2010.mat badlist goodlist gratio
disp('lists of bad/good stations made & saved')

%% Save the M2 amplitudes

count=1;
M2_amp=zeros(length(goodlist),3);
All_p_vals=zeros(length(goodlist),2);
Day_p_vals=zeros(length(goodlist),2);
Night_p_vals=zeros(length(goodlist),2);
tic
for i = 1:length(S);
    load(S(i).name);
    
    % make sure the station had proper data
    if dratio<2/3
        M2_amp(count,1)=All.amplitudes(30);
        M2_amp(count,2)=Day.amplitudes(30);
        M2_amp(count,3)=Night.amplitudes(30);
        
        All_p_vals(count,:)=[All.stats.p(60:61)]';
        Day_p_vals(count,:)=[Day.stats.p(60:61)]';
        Night_p_vals(count,:)=[Night.stats.p(60:61)]';
        
        count=count+1;
    end
end

save 2000-2003_M2_amps_p-vals.mat M2_amp All_p_vals Day_p_vals Night_p_vals
disp('M2 amplitudes extracted & saved')
toc
tic
% Save the O1 amplitudes & p-values

count=1;
O1_amp=zeros(length(goodlist),3);

All_p_vals=zeros(length(goodlist),2);
Day_p_vals=zeros(length(goodlist),2);
Night_p_vals=zeros(length(goodlist),2);

for i = 1:length(S);
    load(S(i).name);
    
    % make sure the station had proper data
    if dratio<2/3
 
        O1_amp(count,1)=All.amplitudes(41);
        O1_amp(count,2)=Day.amplitudes(41);
        O1_amp(count,3)=Night.amplitudes(41);
        
        All_p_vals(count,:)=[All.stats.p(82:83)]';
        Day_p_vals(count,:)=[Day.stats.p(82:83)]';
        Night_p_vals(count,:)=[Night.stats.p(82:83)]';
        
        count=count+1;
    end
end

save 2000-2003_O1_amps_p-vals.mat O1_amp All_p_vals Day_p_vals Night_p_vals
disp('O1 amplitudes extracted & saved')
toc