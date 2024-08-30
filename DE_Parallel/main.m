%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Title: Differential Evolution with Parallel Computing
%%%%%%%%%%%% Writer: Hyeongmeen Baik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clearance and Path
restoredefaultpath

clc;
clear;
close all;
addpath('Parallel_Computing/')
addpath('Differential_Evolution/')

mdl = 'buck_converter_DEMO';
isModelOpen = bdIsLoaded(mdl);


%% Initialization
global Fs Fsamp Ts Tsamp Vg Rdson Tend Ttrig L Cdc Ro


NP =  32;    %Population size: 10N (N - number of dimension)

Kp = 1;
Ki = 1;

CR = 0.9;   %Crossover probability
F = 0.2;    %Differential Weight
ite = 2;

Fs = 10e3;
Fsamp = 10e5;
Ts = 1/Fs;
Tsamp = 1/Fsamp;
Rdson = 0.01;
Tend = 0.1;
Vg = 100;
Ro = 1;
Ttrig = Tend/2;

Vref1 = 40;
Vref2 = 60;

L = 1e-4;
Cdc = 1e-4;

Kp_init = 1;
Ki_init = 1;

% Kp
xi1_min = 0;
xi1_max = 1;
% Ki
xi2_min = 0;
xi2_max = 1;

x_input.min = [xi1_min,xi2_min];
x_input.max = [xi1_max,xi2_max];



%%
DE_out = DE(x_input,NP,CR,F,ite,@MO_DEMO_BUCK_Par);
% x,NP,CR,F,ite,MO

% Termination
if(~isModelOpen)
    close_system(mdl, 0);
end
delete(gcp('nocreate'));


%%
Kp_initial = DE_out.population(1,1,1);
Ki_initial = DE_out.population(1,2,1);
Kp = Kp_initial;
Ki = Ki_initial;
initial_sim = sim('buck_converter_DE_Parallel.slx');
initial_temp.vo = initial_sim.logsout.get('vo').Values;
initial_temp.vref = initial_sim.logsout.get('vref').Values;
initial.time = initial_temp.vo.Time';
initial.vo = initial_temp.vo.Data(:,1)';
initial.vref = initial_temp.vref.Data(:,1)';


Kp_best = DE_out.best_sol(1,1);
Ki_best = DE_out.best_sol(1,1);
Kp = Kp_best;
Ki = Ki_best;
best_sim = sim('buck_converter_DE_Parallel.slx');
best_temp.vo = best_sim.logsout.get('vo').Values;
best_temp.vref = best_sim.logsout.get('vref').Values;
best.time = best_temp.vo.Time';
best.vo = best_temp.vo.Data(:,1)';
best.vref = best_temp.vref.Data(:,1)';

figure(1)
plot(best.time,best.vo,'DisplayName','best');hold on;
plot(initial.time,initial.vo,'DisplayName','initial');
% %% Test
% 
% x = x_input;
% 
% 
% out.N = length(x.max);
% for i = 1:1:NP
%         out.population(i,:,1) = rand(1,out.N).*(x.max-x.min)+x.min;
%         out.y(i,:,1) = MO(out.population(i,:,1));
% end
% k=1;
% for i = 1:1:NP
%     temp1 = randi([1 NP]);
%     temp2 = randi([1 NP]);
%     temp3 = randi([1 NP]);
% 
%     out.donor(i,:) = out.population(temp1,:,k) + F*(out.population(temp2,:,k)-out.population(temp3,:,k));
% 
%     out.CR_comparison_temp = rand([1, out.N])<CR;
%     out.trial(i,:) = out.donor(i,:).*out.CR_comparison_temp +out.population(i,:,k).*~out.CR_comparison_temp;
% 
% end
% 
% x = out.trial;
% 
%     mdl = 'buck_converter_DE_Parallel';
%     isModelOpen = bdIsLoaded(mdl);
%     open_system(mdl);
% 
% 
%     numSims = length(x(:,1));
% 
% for i = 1:1:numSims
%     in_temp(i) = Simulink.SimulationInput(mdl);
%     in(i) = setVariable(in_temp(i),'Kp',x(i,1));
%     in(i) = setVariable(in(i),'Ki',x(i,2));
% end
% 
%     out = parsim(in, 'ShowProgress','on','TransferBaseWorkspaceVariables','on')
% 
% %%
% for i = numSims:-1:1
%         simOut = out(i);
%         vo = simOut.logsout.get('vo').Values;
%         vref = simOut.logsout.get('vref').Values;
%         f(i)=sum((vo.Data(1,:)-vref.Data(1,:)).^2)/length(vo.Data(1,:));
% end
