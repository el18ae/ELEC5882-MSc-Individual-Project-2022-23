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
% The function `FitnessSum` calculates the total fitness value for the given system
% by summing up the individual fitness values for VC1 and IBus. 
%
% This function is used in DPP_PSO Script.
%
%===============================================================================%

function [fitness_result] = FitnessSum(VC1_fitness, w_IBus, IBus_fitness)
    % w_IBus is an unused input after updating the function
    fitness_result = IBus_fitness + VC1_fitness;
end

