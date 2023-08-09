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
% This script models and simulates a PV (Photovoltaic) integrated system
% implementing Differential Power Processing (DPP) from the built simulink model 
% [MUTfinal]. The script is used to analyze and optimize the system using Particle 
% Swarm Optimization (PSO) for selecting optimal values of inductors and capacitors 
% in the circuit.
%
% The model defines parameters for PV irradiance, switching frequency, bus voltage,
% and other related components, and then uses PSO to search the optimal values
% within specified bounds. The simulation is carried out in Simulink, with 
% performance evaluated in terms of fitness based on the waveform analysis method
% on PV1 voltage and IBus current. The results are plotted for visual assessment, 
% providing insights into the system's performance before and after optimization, 
% and the component values are displayed onto the Command Window on Matlab.
%
%===============================================================================%

clear all; clc; close all;
% initialize completion percentage
completion_status = 0;

% Define PV Integrated DPP Simulink Model
G1 = 1000;              % PV1 irradiance
G2 = 1000;              % PV2 irradiance
F = 20e3;               % Switching frequency 
VBus = 36;              % Bus voltage
initial_step = 0.45;    % Simulink "step" block initial value
final_step = 0.5;       % Simulink "step" block final value
step_time = 0.02;       % Simulink "step" block time of step change
tFinal = 0.04;          % Simulink model simulation final time stamp
sample_time = 1e-6;     % 50 samples per period

% weights
w_IBus = 1;

tic
%% PSO
% Define the PSO parameters
n_particles = 16;
max_iter = 22;
w_max = 1; % inertia weight - controls the speed of particle
c1 = 1; % cognitive parameter controls the focus of thier own journey
c2 = 1; % social parameter controls the focus of the collective journey

% Define the search space (bounds) for [L, C]
lb = [1e-6, 1e-6, 1e-3, 1e-3, (1e-6)/2]; % lower bounds
ub = [200e-6, 200e-6, 15e-3, 15e-3, (200e-6)/2]; % upper bounds

% storage initialization
global_best_fit_storage = [];   % global best fitness
global_best_L_storage = [];     % global best L storage
global_best_C_storage = [];     % global best C storage
global_best_Cn_storage = [];     % global best Cn storage
global_best_ripple_storage = [];  % global best ripple
c1_rand_storage = [];
c2_rand_storage = [];
w_rand_storage = [];

% Initialize the particle swarm by creating a struct to store the relevant
% data points. Each struct row will
particle.position = [];
particle.velocity = [];
particle.best.position = [];
particle.best.fitness = [];
particle.best.ripple = [];

% make the struct handle storing data of all particles (n_particles x
% particle) struct:
% B = repmat(A,[M,N]) creates a large matrix B consisting of an M-by-N
% tiling of copies of A, the size of B is[size(A,1)*M, size(A,2)*N].
particles = repmat(particle, n_particles, 1);
global_best.fitness = inf;
global_best.ripple = inf;

% Initialize the particles' positions and velocities
for i = 1:n_particles
    % Generate values from the uniform distribution on the interval (a, b).
    % r = a + (b-a).*rand(row,column);
    particles(i).position = lb + (ub - lb) .* rand(1, 5);   % random position for each particle
    particles(i).velocity = zeros(1, 5);                    % initialize zero velocity for each particle
    C1=particles(i).position(2);
    C2=particles(i).position(1);
    L1=particles(i).position(3);
    L2=particles(i).position(4);
    Cn=particles(i).position(5);
    particles(i).best.position = particles(i).position;     % store the component values

    % store first iteration component values
    C1_OLD = C1;
    C2_OLD = C2;
    L1_OLD = L1;
    L2_OLD = L2;
    Cn_OLD = Cn;

    % Run the Simulink model using the 'sim' command
    sim('MUT_final.slx');

    % extracting relevant values through component
    t = ans.IBus.Time;  % simulation time array
    dt = t(1)-t(2);  % time between each time step in simulation
    VC1 = ans.VC1.Data(:,1);
    VC2 = ans.VC2.Data(:,1);
    IL1 = ans.IL1.Data(:,1);
    IL2 = ans.IL2.Data(:,1);
    VCn = ans.VCn.Data(:,1);
    IBus = ans.IBus.Data(:,1);

    % compute fitness
    [VC1_fitness,~] = FitnessRunVC1(C1_,C2_,step_time,t,VC1);
    [IBus_fitness] = FitnessRunIBus(C1_,C2_,step_time,t,IBus);

    % store the fitness
    [particles(i).best.fitness] = FitnessSum(VC1_fitness, w_IBus, IBus_fitness);
    fitness_OLD = particles(i).best.fitness;

    % check which particle had the best fitness compared to all (collective journey)
    if particles(i).best.fitness < global_best.fitness
        global_best = particles(i).best;
    end
end

% Perform the PSO iterations
for iter = 1:max_iter
    w = w_max;
    for i = 1:n_particles
        % create randomly generated inertia weight, cognetive and social
        % parameters
        w_rand = w * (0.01+(1-0.01).*rand(1,1));  % random value between 0.01 & 1
        c1_rand = c1 * rand(1, 5);  % 1x5 matrix of random values between 0 & 1 
        c2_rand = c2 * rand(1, 5);

        % store them
        c1_rand_storage = [c1_rand_storage w_rand];
        c2_rand_storage = [c2_rand_storage c1_rand];
        w_rand_storage = [w_rand_storage c2_rand];

        % Update the particle's velocity
        particles(i).velocity = w_rand * particles(i).velocity ...
            + c1_rand .* (particles(i).best.position - particles(i).position) ...
            + c2_rand .* (global_best.position - particles(i).position);

        % Update the particle's position
        particles(i).position = particles(i).position + particles(i).velocity;

        % Apply bounds - do not exceed
        particles(i).position = max(particles(i).position, lb);
        particles(i).position = min(particles(i).position, ub);

        % Evaluate the fitness of the new position
        C1 = particles(i).position(2);
        C2 = particles(i).position(1);
        L1 = particles(i).position(3);
        L2 = particles(i).position(4);
        Cn = particles(i).position(5);

        % Run the Simulink model using the 'sim' command
        sim('MUT_final.slx');

        % completion status display
        clc;
        completion_status = completion_status + 1;
        status_percent = (completion_status/(max_iter*n_particles))*100;
        disp(['completion = ' num2str(status_percent) '%']);

        % extracting relevant values through component
        t = ans.IBus.Time;  % simulation time array
        dt = t(1)-t(2);  % time between each time step in simulation
        VC1 = ans.VC1.Data(:,1);
        VC2 = ans.VC2.Data(:,1);
        IL1 = ans.IL1.Data(:,1);
        IL2 = ans.IL2.Data(:,1);
        VCn = ans.VCn.Data(:,1);
        IBus = ans.IBus.Data(:,1);

        % compute fitness
        [VC1_fitness,~] = FitnessRunVC1(C1_,C2_,step_time,t,VC1);
        [IBus_fitness] = FitnessRunIBus(C1_,C2_,step_time,t,IBus);

        % store the fitness
        [particles(i).fitness] = FitnessSum(VC1_fitness, w_IBus, IBus_fitness);

        % Update the particle's personal best (own journey)
        if particles(i).fitness < particles(i).best.fitness
            particles(i).best.position = particles(i).position;
            particles(i).best.fitness = particles(i).fitness;

            % Update the global best  (collective journey)
            if particles(i).best.fitness < global_best.fitness
                global_best = particles(i).best;
            end
        end
    end

    % store results
    global_best_fit_storage = [global_best_fit_storage global_best.fitness];
    global_best_L_storage = [global_best_L_storage global_best.position(3)];
    global_best_C_storage = [global_best_C_storage global_best.position(2)];
    global_best_Cn_storage = [global_best_Cn_storage global_best.position(5)];
end
toc
% extract minimum fitness
fitness_best = min(global_best_fit_storage);
index = find(global_best_fit_storage == fitness_best, 1);

% Print the final results [C1,C2,Cn,L1,L2,old_VC1_fitness,old_VC1_ripple];
disp('AFTER OPTIMISATION');
disp(['Minimum Fitness = ' num2str(min(global_best_fit_storage))]);
disp(['Optimal L1 = ' num2str(global_best.position(3))]);
disp(['Optimal L2 = ' num2str(global_best.position(4))]);
disp(['Optimal C1 = ' num2str(global_best.position(2))]);
disp(['Optimal C2 = ' num2str(global_best.position(1))]);
disp(['Optimal Cn = ' num2str(global_best.position(5))]);

% initialize component values with optimal values
C1 = global_best.position(2);
C2 = global_best.position(1);
L1 = global_best.position(3);
L2 = global_best.position(4);
Cn = global_best.position(5);

% Run the Simulink model using the 'sim' command
sim('MUT_final.slx');

% extracting simulation relevant values through component
t = ans.IBus.Time;  % simulation time array
dt = t(1)-t(2);  % time between each time step in simulation
VC1 = ans.VC1.Data(:,1);
VC2 = ans.VC2.Data(:,1);
IL1 = ans.IL1.Data(:,1);
IL2 = ans.IL2.Data(:,1);
VCn = ans.VCn.Data(:,1);
IBus = ans.IBus.Data(:,1);

% plot after optimisation
figure
plot(t, VC1,'-','color','k','LineWidth',2,'DisplayName',['VC1']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VC1 - after optimisation']);
hold on;
grid on;

figure
plot(t, VC2,'-','color','k','LineWidth',2,'DisplayName',['VC2']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VC2 - after optimisation']);
hold on;
grid on;

figure
plot(t, IL1,'-','color','k','LineWidth',2,'DisplayName',['IL1']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['IL1 - after optimisation']);
hold on;
grid on;

figure
plot(t, IL2,'-','color','k','LineWidth',2,'DisplayName',['IL2']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['IL2 - after optimisation']);
hold on;
grid on;

figure
plot(t, VCn,'-','color','k','LineWidth',2,'DisplayName',['VCn']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VCn - after optimisation']);
hold on;
grid on;

figure
plot(t, IBus,'-','color','k','LineWidth',2,'DisplayName',['IBus']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['IBus - after optimisation']);
hold on;
grid on;

% plot change of fitness
figure
plot(1:max_iter, global_best_fit_storage);
title('Change in global best fitness over iterations');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation of first iteration random component values 
C1 = C1_OLD;
C2 = C2_OLD;
L1 = L1_OLD;
L2 = L2_OLD;
Cn = Cn_OLD;

disp('BEFORE OPTIMISATION');
disp(['Minimum Fitness = ' num2str(fitness_OLD)]);
disp(['Optimal L1 = ' num2str(L1_OLD)]);
disp(['Optimal L2 = ' num2str(L2_OLD)]);
disp(['Optimal C1 = ' num2str(C1_OLD)]);
disp(['Optimal C2 = ' num2str(C2_OLD)]);
disp(['Optimal Cn = ' num2str(Cn_OLD)]);

% Run the Simulink model using the 'sim' command
sim('MUT_final.slx');

% extracting relevant values through component
t = ans.IBus.Time;  % simulation time array
dt = t(1)-t(2);  % time between each time step in simulation
VC1 = ans.VC1.Data(:,1);
VC2 = ans.VC2.Data(:,1);
IL1 = ans.IL1.Data(:,1);
IL2 = ans.IL2.Data(:,1);
VCn = ans.VCn.Data(:,1);
IBus = ans.IBus.Data(:,1);

% plot before optimisation
figure
plot(t, VC1,'-','color','k','LineWidth',2,'DisplayName',['VC1']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VC1 - before optimisation']);
hold on;
grid on;

figure
plot(t, VC2,'-','color','k','LineWidth',2,'DisplayName',['VC2']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VC2 - before optimisation']);
hold on;
grid on;

figure
plot(t, IL1,'-','color','k','LineWidth',2,'DisplayName',['IL1']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['IL1 - before optimisation']);
hold on;
grid on;

figure
plot(t, IL2,'-','color','k','LineWidth',2,'DisplayName',['IL2']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['IL2 - before optimisation']);
hold on;
grid on;

figure
plot(t, VCn,'-','color','k','LineWidth',2,'DisplayName',['VCn']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['VCn - before optimisation']);
hold on;
grid on;

figure
plot(t, IBus,'-','color','k','LineWidth',2,'DisplayName',['IBus']);
xlabel('TIME (s)');
ylabel('Voltage (V)');
title(['IBus - before optimisation']);
hold on;
grid on;
