function [h]=F_obs_psd(data,t, startd, endd,name)
close all;
ind1=find(t==startd);
ind2=find(t==endd);

t= t(ind1:ind2);
data= data(ind1:ind2);

test=isnan(data);
dratio=sum(test)/length(data)

% if the ratio is good

if dratio < 0.1
    
    dataI=interp1(t,data,t,'cubic');
    
    NW=4; %2, 5/2, 3, 7/2, 4 are the typical values
    NFFT=6*600*60; % 1 month window
    Fr=1/60; % Sampling frequency
    % Determining the PSD
    [PX,PXc,fx]=pmtm(dataI,NW,NFFT,Fr);
    % Determining the power amplitudes for Xd
    AmpX=zeros(length(PX),1);
    for i=1:length(PX)
        AmpX(i,1)=(2*PX(i,1)*(Fr/NFFT));
    end
    % Determining the periods
    TX=zeros(length(fx)-1,1);
    for i=2:length(fx)
        TX(i-1)=(1/fx(i))*(1/60)*(1/60);
    end

    % for plotting
    periods=[4 4.8 6 8 11.967236 12 12.42060122 12.6583 ... 
         23.934472 24 24.0659 25.8913 25.8956];
    period_string=['S6';'S5';'S4';'S3';'K2';'S2';'M2';'N2';...
         'K1';'S1';'P1';'O1';'Q1'];
     
    h=figure(2);
    semilogy(TX,PX(2:end),'k','LineWidth',1)
    set(gca,'FontSize',16)
    xlim([0 30])

    % Plotting the tidal lines
    a = axis;
    for kk=1:length(periods)
        line([periods(kk) periods(kk)],[a(3) a(4)],'LineStyle','--','color','b') 
        % Plotting their labels
        text(periods(kk),kk^7-100,period_string(kk,1:2))
    end
    xlabel('Period (Hours)')
    ylabel('Power/Frequency (nT^2/Hz)')

    year1=datestr(startd);
    year1=year1(8:11);
    year2=datestr(endd);
    year2=year2(8:11);

    figname=[name '_psd_' year1 '-' year2 '.eps']
    print(figname,'-deps')
    
elseif dratio >= 0.1
    h= NaN
    disp('nope')
end