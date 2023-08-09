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
% This function computes the fitness and ripple for VC1 based on the input parameters, 
% determining the reference value for the current, computing the time when the current 
% reaches rise time final point, and calculating the areas above and below the reference 
% line. A weighted sum of these areas is then used to compute the fitness, considering 
% transient and steady-state behavior, and the ripple is also determined.
%
% This function is used in DPP_PSO Script.
%
%===============================================================================%

function [fitness_result, ripple_result] = FitnessRunVC1(C1_,C2_,step_time,t,IL)
% Fitness and Ripple storage
fitness = zeros(length(C1_));
rippleIL = zeros(length(C1_));

%% Loop
for i = 1:length(C1_)

% Find the location of the last occurrence of the value 0.05 in the time data
idx_last_occurrence = max(find(t == step_time));

% Extract the time and current data up to the last occurrence of 0.02
t_subset = t(1:idx_last_occurrence);
IL_subset = IL(1:idx_last_occurrence);

% Define the fraction of the subset to use for the reference line
frac_ref = 0.1;

% Compute the reference value as the average of the last frac_ref of the subset
ref_value = mean(IL_subset(end-round(frac_ref*numel(IL_subset))+1:end));

% Compute the time when the current reaches 90% of its maximum value
IL_max = ref_value;
flag1 = 0;  % flag to check first Bref crossing

% check Bref crossing to stop right before the steady state
for cc = 1:length(IL_subset)
    if (IL_subset(cc) >= IL_max && flag1 == 0)
        idx_90 = cc;
        flag1 = 1;
    end
end

% Extract the time and current data up to the risetime final point
IL_rise = IL(1:idx_90);
t_rise = t(1:length(IL_rise));

% Compute the area above and below the reference line for the rise data
% argument (0) sets negative values of computed area to zero as we only
% interested in max alone and min alone.
area_rise = abs((IL_rise - ref_value));

% Compute the area above and below the reference line for the remaining data
t_rest = t(length(t_rise)+1:end);
IL_rest = IL(length(t_rise)+1:end);
area_rest = abs((IL_rest - ref_value));

% Assign a higher weight to the areas computed after the 90% point
%w = 0.8;
w1 = 30;  % should be higher (not too high) than w2 for improvement of transient(focus)
w2 = 1;
% optimisation of optimisation

% Compute the fitness as a weighted sum of the areas then divide by total
% time taken to normalize
fitness(i) = ((w1)*(sum(area_rise)) + (w2)*(sum(area_rest)))/max(t);

% Ripple computation
Ipp = max(IL_subset(end-round(frac_ref*numel(IL_subset))+1:end))-min(IL_subset(end-round(frac_ref*numel(IL_subset))+1:end));  % Ipk-pk
Iavg = ref_value;  % average current
rippleIL(i) = 100*(Ipp/Iavg);
ripple_result = rippleIL(i);
fitness_result = fitness(i);
end
end

