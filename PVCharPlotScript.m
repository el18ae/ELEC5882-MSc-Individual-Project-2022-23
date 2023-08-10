%{
    NAME: ABDULLAH ESSA
    SID: 201256467
    ELEC5564M - Electric Power Generation by Renewable Sources
    Assignment (1), Question (1) - Write a MATLAB program to solve the 
    mathematical model of a solar PV array 

    Last Updated: 25 July 2023
%}
clear; clc; close all;
%% 1.1) At an ambient temperature (Ta) of 25°C
% (i) Plot the current-voltage (I-V) and power-voltage (P-V) 
% characteristics of a solar PV array for different irradiance levels 
% G=0.3, 0.5 and 1.0 kW/m2.

% initializing operational conditions
G = [0.3 0.5 1];  % irradiation
Ta = [25];  % ambient temperature
Rp = 5e+5;
% Simulating & Plotting
FF = PVComp(G, Ta,Rp,1);

%(ii) Calculate the corresponding fill factor for each irradiance level and 
% explain how the characteristics of a solar PV array e.g., (Isc(Isc​, 
% VocVoc​, VMPPVMPP​, IMPPIMPP​ and MPP)MPP) may vary with varying irradiance 
% levels.  

% explanation of how the characteristics of solar PV array may vary is
% found in the submitted PDF
figure(2);
% initialize plot color array
C = {'b','r',"#77AC30","#4DBEEE",'m'};
% plot the FF against irradiance
for i = 1:length(FF)
    scatter(G(i), FF(i),"MarkerFaceColor", C{i},'DisplayName',['(G = ',num2str(G(i)),')&(FF = ',num2str(FF(i)),')']);
    xlabel('Irradiance (G)');
    ylabel('Fill Factor');
    grid on;
    legend('-DynamicLegend','Location','NorthWest');
    hold on;
end
plot(G, FF, 'k','DisplayName',['Line']);
title(['Fill Factor - for Ta(',num2str(Ta),') & G =(', num2str(G),')']);
grid on;

% (iii) answer found in the submitted PDF file 

%% 1.2) At irradiance level of G=1.0 kW/m2
% (i) Plot the current-voltage (I-V) and power-voltage (P-V) 
% characteristics of a solar PV array for different range of 
% temperatures Ta of 25, 40 and 60°C

% initializing operational conditions
G = [1];  % irradiation
Ta = [25 40 60];  % ambient temperature
Rp = 5e+5;
% Simulating & Plotting
FF = PVComp(G, Ta,Rp,3);

%(ii) Calculate the corresponding fill factor on each Temperature level and 
% explain how the characteristics of a solar PV array e.g., (Isc(Isc​, 
% VocVoc​, VMPPVMPP​, IMPPIMPP​ and MPP)MPP) may vary with varying temperature 
% levels.  

% explanation of how the characteristics of solar PV array may vary is
% found in the submitted PDF
figure(4);
% plot the FF against irradiance
for i = 1:length(FF)
    scatter(Ta(i), FF(i),"MarkerFaceColor",C{i},'DisplayName',['(Ta = ',num2str(Ta(i)),')&(FF = ',num2str(FF(i)),')']);
    xlabel('Temperature (Ta)');
    ylabel('Fill Factor');
    grid on;
    legend('-DynamicLegend','Location','NorthWest');
    hold on;
end
plot(Ta, FF, 'k','DisplayName',['Line']);
title(['Fill Factor - for Ta(',num2str(Ta),') & G =(', num2str(G),')']);
grid on;

% (iii) answer found in the submitted PDF file 

%% 1.3) Code submission of aforementioned in submission link on minerva.

%% 1.4) Analyzing the affect of varying Rp on constant G and Ta
% At a constant irradiance level, G of 1.0 kW/m2  and an ambient 
% temperature (Ta) of  25°C,  plot and study the impact of decreasing the 
% value of Rp from 5×10^5 Ω  to  0.5,  0.05  and 0.005Ω on the 
% corresponding MPP, Isc​ and Voc​ and comment on the results. 

% initializing operational conditions
G = [1];  % irradiation
Ta = [25];  % ambient temperature
Rp = [5e+5 0.5 0.05 0.005];

% Simulating & Plotting
FF = PVComp(G, Ta,Rp,5);
