%===============================================================================%
% ELEC5882 MSc Individual Project 2022/23
%===============================================================================%
%
% Name: Abdullah Essa
% Student ID: 201256467
% University: University of Leeds, School of Electrical and Electronics
% Supervisor: Dr. Benjamin Chong
% Last updated: 09 August 2023
%
%===============================================================================%
%
% Description:
%
% This script sets up and runs a simulation for a PV (Photovoltaic) system
% with DPP, defining various component values and parameters. The code is prepared
% to handle both continuous and single step changes in the system. For a continuous
% variation, 'MUTfinal2.slx' is run, and for a single step change, 'MUTfinal.slx'
% can be used by uncommenting the relevant sections.
%
% In the single step change version, it computes the rise/fall time and ripple for
% voltage VC1 and current IBus, displaying them in the command window.
% Plots for voltage (VC1, VC2) and current (IBus) are generated, showing the behavior
% after optimization, with figures labeled and axis limits set for clear representation.
%
%===============================================================================%

clear; clc; close all;

%% DPP Initialisation
% define model values
C1_ = 1.59e-6;
C2_ = 1.59e-6;
C1 = C1_;
C2 = C2_;
Cn = (1.59/2)*1e-6;
L1 = 4.9e-3;
L2 = 4.9e-3;
G1 = 1000;
G2 = 1000;
F = 20e3;
VBus = 36;
initial_step = 0.45;
final_step = 0.5;
tFinal = 0.14;
step_time = 0.02;
sample_time = 1e-6;  % 50 samples per period

% Run the Simulink model for PV DPP with various step changes
sim('MUTfinal2.slx');

% uncomment below if PV DPP Simulink for a single step change
% sim('MUTfinal.slx');

% extracting relevant values through component
t = ans.IBus.Time;  % simulation time array
dt = t(1)-t(2);  % time between each time step in simulation
VC1 = ans.VC1.Data(:,1);
VC2 = ans.VC2.Data(:,1);
IL1 = ans.IL1.Data(:,1);
IL2 = ans.IL2.Data(:,1);
VCn = ans.VCn.Data(:,1);
IBus = ans.IBus.Data(:,1);

% uncomment below if PV DPP Simulink for a single step change
% [ripple_vc1, ripple_ibus,rise_time_vc1,rise_time_ibus] = GetRipple(initial_step,final_step,step_time,t,VC1,IBus);
%
% disp(['Voltage VC1 Rise/Fall Time = ' num2str(rise_time_vc1)]);
% disp(['Voltage VC1 Ripple = ' num2str(ripple_vc1)]);
% disp(['Current IBus Ripple = ' num2str(ripple_ibus)]);
% disp(['Current IBus Rise/Fall Time = ' num2str(rise_time_ibus)]);

figure
plot(t, VC1,'-','color','k','LineWidth',2,'DisplayName',['VC1']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VC1 - after optimisation']);
hold on;
grid on;
axis([0.005 0.14 17.5 20.5]);
% uncomment and edit below if PV DPP Simulink for a single step change
% axis([0.005 0.02 17.5 20.5]);

figure
plot(t, VC2,'-','color','k','LineWidth',2,'DisplayName',['VC2']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VC2 - after optimisation']);
hold on;
grid on;
axis([0.005 0.14 15.5 19]);
% uncomment and edit below if PV DPP Simulink for a single step change
% axis([0.005 0.02 15.5 19]);

figure
plot(t, IBus,'-','color','k','LineWidth',2,'DisplayName',['IBus']);
xlabel('TIME (s)');
ylabel('Current (A)');
title(['IBus - after optimisation']);
hold on;
grid on;
axis([0.005 0.14 2 5.5]);
% uncomment and edit below if PV DPP Simulink for a single step change
% axis([0.005 0.02 2 5.5]);