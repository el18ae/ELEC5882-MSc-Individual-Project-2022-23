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
% This script creates a grouped bar graph that visualizes the comparison
% between two approaches, 'Estimation Equation' and 'Optimization', for
% two different measurements: 'Vpv1 measurement' and 'IBus measurement'.
% The blue bars represent the 'Estimation Equation' method, and the red
% bars represent the 'Optimization' method.
%
%===============================================================================%

clear; clc; close all;

% Define data that will be compared:
% Rows represent measurements ('Vpv1 measurement' and 'IBus measurement')
% Columns represent approaches ('Estimation Equation' and 'Optimization')
data = [1.03 0.025; 3 0.091]; 

% Define categories
categories = {'Vpv1 measurement', 'IBus measurement'}; 

% Create a figure
figure;

% Create bar graph with grouped bars
b = bar(data, 'grouped');

% Set categories
set(gca, 'XTickLabel', categories);

% Set colors for each group
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';

% Add legend
legend('Estimation Equation', 'Optimization');

% Add title
title('Comparison of Ripple (%) Between Estimation Equation and Optimization Approaches');

% Label y axis
ylabel('Ripple (%)');

ylim([0, 3.5]);
