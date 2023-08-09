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
% The function calculates a fitness value for the IBus based on the given parameters. 
% The calculation is divided into two parts, where each part involves the extraction 
% of subsets of the time and current data and the computation of fitness values using 
% weighted areas above and below a reference line.
%
% This function is used in DPP_PSO Script.
%
%===============================================================================%

function [fitness_result] = FitnessRunIBus(C1_,C2_,step_time,t,IL)
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

% check Bref crossing to stop right before the steady state rise time point
for cc = 2:length(IL_subset)
    if (IL_subset(cc) <= IL_max && flag1 == 0)
        idx_90 = cc;
        flag1 = 1;
    end
end

% Extract the time and current data up to the rise time point
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

% Assign a higher weight to the areas computed after the rise time point
%w = 0.8;
w1 = 30;  % should be higher (not too high) than w2 for improvement of transient(focus)
w2 = 1;
% optimisation of optimisation

% Compute the fitness as a weighted sum of the areas then divide by total
% time taken to normalize
fitness(i) = ((w1)*(sum(area_rise)) + (w2)*(sum(area_rest)))/max(t);
fitness_result1 = fitness(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the location of the last occurrence of the value 0.05 in the time data
idx_final_occurrence = max(find(t == max(t)));

% Extract the time and current data up to the last occurrence of 0.02
t_final_subset = t(idx_last_occurrence:idx_final_occurrence);
IL_final_subset = IL(idx_last_occurrence:idx_final_occurrence);

% Define the fraction of the subset to use for the reference line
frac_ref = 0.1;

% Compute the reference value as the average of the last frac_ref of the subset
ref_value_final = mean(IL_final_subset(end-round(frac_ref*numel(IL_final_subset))+1:end));

% Compute the time when the current reaches rise time point
IL_max_final = ref_value_final;
flag1_final = 0;  % flag to check first Bref crossing

% check Bref crossing to stop right before the steady state
for cc = 1:length(IL_final_subset)
    if (IL_final_subset(cc) >= IL_max_final && flag1_final == 0)
        idx_90_final = cc;
        flag1_final = 1;
    end
end

% Extract the time and current data up to the rise time point
IL_rise_final = IL(idx_last_occurrence:idx_90_final);
t_rise_final = t(idx_last_occurrence:length(IL_rise_final));

% Compute the area above and below the reference line for the rise data
% argument (0) sets negative values of computed area to zero as we only
% interested in max alone and min alone.
area_rise_final = abs((IL_rise_final - ref_value_final));

% Compute the area above and below the reference line for the remaining data
t_rest_final = t(length(t_rise_final)+1:end);
IL_rest_final = IL(length(t_rise_final)+1:end);
area_rest_final = abs((IL_rest_final - ref_value_final));

% Assign a higher weight to the areas computed after the 90% point
%w = 0.8;
w1 = 30;  % should be higher (not too high) than w2 for improvement of transient(focus)
w2 = 1;
% optimisation of optimisation

% Compute the fitness as a weighted sum of the areas then divide by total
% time taken to normalize
fitness(i) = ((w1)*(sum(area_rise_final)) + (w2)*(sum(area_rest_final)))/max(t);
fitness_result2 = fitness(i);
fitness_result = fitness_result1 + fitness_result2;
end
end

