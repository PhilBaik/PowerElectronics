%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Title: Parallel Computing (DEMO2)
%%%%%%%%%%%% Writer: Hyeongmeen Baik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% For buck converter, with block parameters

clc;
clear;
close all;

Fswit = 10e3;
Tswit = 1/Fswit;
Tend = 10;
R = 10;
D = 0.5;

mdl = 'buck_converter_DEMO';
isModelOpen = bdIsLoaded(mdl);
open_system(mdl);

R_sweep = R*(1:1:10);
numSims = length(R_sweep);

%% by For loop

start_single = tic;
for i = 1:1:numSims
R_single = R_sweep(1,i);
in_temp_single = Simulink.SimulationInput(mdl);
in_single = setBlockParameter(in_temp_single,[mdl '/Buck_converter'], 'R', num2str(R_single));
out_single(i) = sim(in_single);
disp_name = sprintf('iteration: %2.f %%',i/numSims*100)
end
T_single = toc(start_single);
%% by Parallel computing
start_par = tic;
for i = 1:1:numSims
    in_temp(i) = Simulink.SimulationInput(mdl);
    in(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'R', num2str(R_sweep(i)));
end
out = parsim(in, 'ShowProgress','on','TransferBaseWorkspaceVariables','on');
T_par = toc(start_par);

%%
k_improve = T_single/T_par;
%%
% Termination
if(~isModelOpen)
    close_system(mdl, 0);
end
delete(gcp('nocreate'));
%% Plot
figure(11)
for i = numSims:-1:1
        simOut = out(i);
        ts = simOut.logsout.get('vo').Values;
        % 'ts' is a MATLAB 'timeseries' object that stores the time and
        % data values for the logged 'vertical_disp' signal.
        % Use the 'plot' method of the object to plot the data against the
        % time.
        plot(ts);
        legend_labels{i} = ['Run ' num2str(i)];
        hold all
end
title('Voltage Output')
xlabel('Time (s)');
ylabel('Voltage (V)');
legend(legend_labels,'Location','NorthEastOutside');

figure(12)
for i = numSims:-1:1
        simOut = out(i);
        ts = simOut.logsout.get('io').Values;
        % 'ts' is a MATLAB 'timeseries' object that stores the time and
        % data values for the logged 'vertical_disp' signal.
        % Use the 'plot' method of the object to plot the data against the
        % time.
        plot(ts);
        legend_labels{i} = ['Run ' num2str(i)];
        hold all
end
title('Current Output')
xlabel('Time (s)');
ylabel('Current (A)');
legend(legend_labels,'Location','NorthEastOutside');


%% For buck converter, with Simulation Input


Fswit = 10e3;
Tswit = 1/Fswit;
Tend = 1;
R = 10;
D = 0.5;

mdl = 'buck_converter_DEMO';
isModelOpen = bdIsLoaded(mdl);
open_system(mdl);


D_sweep = 0.1:0.1:0.9; 
numSims = length(D_sweep);

%%
for i = 1:1:numSims
    in_temp(i) = Simulink.SimulationInput(mdl);
    in(i) = setVariable(in_temp(i),'D', D_sweep(i));
end
out = parsim(in, 'ShowProgress','on','TransferBaseWorkspaceVariables','on');


%% Plot
figure(11)
for i = numSims:-1:1
        simOut = out(i);
        ts = simOut.logsout.get('vo').Values;
        % 'ts' is a MATLAB 'timeseries' object that stores the time and
        % data values for the logged 'vertical_disp' signal.
        % Use the 'plot' method of the object to plot the data against the
        % time.
        plot(ts);
        legend_labels{i} = ['Run ' num2str(i)];
        hold all
end
title('Voltage Output')
xlabel('Time (s)');
ylabel('Voltage (V)');
legend(legend_labels,'Location','NorthEastOutside');

figure(12)
for i = numSims:-1:1
        simOut = out(i);
        ts = simOut.logsout.get('io').Values;
        % 'ts' is a MATLAB 'timeseries' object that stores the time and
        % data values for the logged 'vertical_disp' signal.
        % Use the 'plot' method of the object to plot the data against the
        % time.
        plot(ts);
        legend_labels{i} = ['Run ' num2str(i)];
        hold all
end
title('Current Output')
xlabel('Time (s)');
ylabel('Current (A)');
legend(legend_labels,'Location','NorthEastOutside');


