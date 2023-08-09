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
% This function computes the rise time and ripple of the IBus current and VC1 
% voltage via waveform analysis. It calculates the ripple and rise time under 
% different scenarios (increasing/decreasing duty ratio or no change) and returns 
% the values for both IBus and VC1.
%
%===============================================================================%

function [ripple_vc1, ripple_ibus,rise_time_vc1,rise_time_ibus] = GetRipple(initial_step,final_step,step_time,t,VC1,IBus)

% Find the location of the last occurrence of the value 0.02 in the time data
idx_last_occurrence = max(find(t == step_time));

% Extract the time and current data from 0 to the last occurrence of 0.02
t_subset = t(1:idx_last_occurrence);

% Find the location of the last occurrence of the value 0.04 in the time data
idx_final_occurrence = max(find(t == max(t)));

% Extract the time and current data fromm 0.02 to the last occurrence of 0.04
t_final_subset = t(idx_last_occurrence:idx_final_occurrence);

% Define the fraction of the subset to use for the reference line
frac_ref = 0.1;

%% IBus ripple
IBus_first_subset = IBus(1:idx_last_occurrence);
IBus_second_subset = IBus(idx_last_occurrence:idx_final_occurrence);

% Compute the reference value as the average of the last frac_ref of the subset
IBus_secondsub_ref_value = mean(IBus_second_subset(end-round(frac_ref*numel(IBus_second_subset))+1:end));

% Ripple computation
Ipp = max(IBus_second_subset(end-round(frac_ref*numel(IBus_second_subset))+1:end))-min(IBus_second_subset(end-round(frac_ref*numel(IBus_second_subset))+1:end));  % Ipk-pk
Iavg = IBus_secondsub_ref_value;  % average current
ripple_ibus = 100*(Ipp/Iavg);

disp(['Current IBus Average = ' num2str(Iavg)]);

% rise time computation
flagi_final = 0;  % flag to check first Bref crossing

% rise time computation
if(initial_step < final_step)  % if K_i < K_i+1 (increasing duty ratio)

    % check Bref crossing to stop right before the steady state
    for cc = 1:length(IBus_second_subset)
        if (IBus_second_subset(cc) >= IBus_secondsub_ref_value && flagi_final == 0)
            idx_RISEi = cc;
            flagi_final = 1;
            break;
        end
    end
elseif(initial_step > final_step)  % if K_i > K_i+1 (decreasing duty ratio)

    for cc = 1:length(IBus_second_subset)
        if (IBus_second_subset(cc) <= IBus_secondsub_ref_value && flagi_final == 0)
            idx_RISEi = cc;
            flagi_final = 1;
            break;
        end
    end

else  % when both steps are similar and no change
    flagi_final = 0;  % set it back to zero to indicate no rise time detection
end

if (flagi_final == 1)
    rise_time_ibus = t_final_subset(idx_RISEi)-t_final_subset(1);
else
    rise_time_ibus = 0;
end

%% VC1 ripple
VC1_first_subset = VC1(1:idx_last_occurrence);
VC1_second_subset = VC1(idx_last_occurrence:idx_final_occurrence);

% Compute the reference value as the average of the last frac_ref of the subset
VC1_secondsub_ref_value = mean(VC1_second_subset(end-round(frac_ref*numel(VC1_second_subset))+1:end));

% Ripple computation
Vpp = max(IBus_second_subset(end-round(frac_ref*numel(IBus_second_subset))+1:end))-min(IBus_second_subset(end-round(frac_ref*numel(IBus_second_subset))+1:end));  % Vpk-pk
Vavg = VC1_secondsub_ref_value;  % average voltage
ripple_vc1 = 100*(Vpp/Vavg);
disp(['Voltage VC1 Average = ' num2str(Vavg)]);

flagV_final = 0;  % flag to check first Bref crossing

% rise time computation
if(initial_step < final_step)  % if K_i < K_i+1 (increasing duty ratio)

    % check Bref crossing to stop right before the steady state
    for cc = 1:length(VC1_second_subset)
        if (VC1_second_subset(cc) <= VC1_secondsub_ref_value && flagV_final == 0)
            idx_RISEV = cc;
            flagV_final = 1;
            break;
        end
    end
elseif(initial_step > final_step)  % if K_i > K_i+1 (decreasing duty ratio)

    for cc = 1:length(VC1_second_subset)
        if (VC1_second_subset(cc) >= VC1_secondsub_ref_value && flagV_final == 0)
            idx_RISEV = cc;
            flagV_final = 1;
            break;
        end
    end

else  % when both steps are similar and no change
    flagV_final = 0;  % set it back to zero to indicate no rise time detection
end

if (flagV_final == 1)
    rise_time_vc1 = t_final_subset(idx_RISEV)-t_final_subset(1);
else
    rise_time_vc1 = 0;
end

end

