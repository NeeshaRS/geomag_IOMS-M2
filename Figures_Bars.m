clc; clear all; close all;

load stations_used_list_1995-1998.mat
% load M2_amps_p-vals.mat
load O1_amps_p-vals.mat
%% M2 bar plots
n=length(goodlist);

% for i=1
for i=1:n
    figure(1); 
    hold on;
    set(gca, 'FontSize', 21, 'LineWidth', 2); 
    y=M2_amp(i,:);
    h1= bar(1,y(1)); col = 'k'; set(h1, 'FaceColor', col, 'LineWidth', 2); 
    h2= bar(2,y(2)); col = 'r'; set(h2, 'FaceColor', col, 'LineWidth', 2);
    h3= bar(3,y(3)); col = 'b'; set(h3, 'FaceColor', col, 'LineWidth', 2);
    ylim([0 6.5])
    set(gca,'XTickLabel',{'All','Day','Night'},'XTick',[1 2 3 ])
    ylabel('Amplitude (nT)')
    
    picname=['1995-1998_Solar_Minimum/M2_Bar_Plots/Same_Scale/' goodlist(i,:) '_bar_plot.png'];
    print(picname,'-dpng');
    close all
end

disp('M2 bar plots made')

%% O1 bar plots
n=length(goodlist);
tic
% for i=1
for i=1:n
    figure(1); 
    hold on;
    set(gca, 'FontSize', 21, 'LineWidth', 2); 
    y=O1_amp(i,:);
    h1= bar(1,y(1)); col = 'k'; set(h1, 'FaceColor', col, 'LineWidth', 2); 
    h2= bar(2,y(2)); col = 'r'; set(h2, 'FaceColor', col, 'LineWidth', 2);
    h3= bar(3,y(3)); col = 'b'; set(h3, 'FaceColor', col, 'LineWidth', 2);
%     ylim([0 4.6])
    set(gca,'XTickLabel',{'All','Day','Night'},'XTick',[1 2 3 ])
    ylabel('Amplitude (nT)')
    
    picname=['1995-1998_Solar_Minimum/O1_Bar_Plots/Individual_Scale/' goodlist(i,:) '_bar_plot.png'];
    print(picname,'-dpng');
    close all
end

disp('O1 bar plots made')
toc