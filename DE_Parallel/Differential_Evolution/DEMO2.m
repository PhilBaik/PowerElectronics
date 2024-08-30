%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Title: Differential Evolution (DEMO2)
%%%%%%%%%%%% Writer: Hyeongmeen Baik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% Initialization

global Fs Fsamp Ts Tsamp Vg Rdson Tend Kp Ki Ttrig L Cdc Ro

xi1_min = 0;
xi1_max = 1;
xi2_min = 0;
xi2_max = 1;

NP = 20 ;    %Population size: 10N (N - number of dimension)

x_input.min = [xi1_min,xi2_min];
x_input.max = [xi1_max,xi2_max];

CR = 0.9;   %Crossover probability
F = 0.2;    %Differential Weight
ite = 10;

Fs = 10e3;
Fsamp = 10e5;
Ts = 1/Fs;
Tsamp = 1/Fsamp;
Rdson = 0.01;
Tend = Ts*1000;
Vg = 100;
Ro = 1;

Vref1 = 40;
Vref2 = 60;

L = 1e-4;
Cdc = 1e-4;

Kp_init = 1;
Ki_init = 1;


%% Simulink Test
n_cycles = 500;
Tstart = Tend - Ts*n_cycles;
Ttrig = (Tstart+Tend)/2;
n_points = round(n_cycles*Ts/Tsamp);

Kp = Kp_init;
Ki = Ki_init;
temp_sim_out = sim('buck_demo.slx');

sim_out.time = temp_sim_out.sim_out.Time';
sim_out.time(1,end-n_points)
sim_out.vo = temp_sim_out.sim_out.Data(:,1)';
sim_out.vref = temp_sim_out.sim_out.Data(:,2)';

sim_out.time = sim_out.time(1,end-n_points:end);
sim_out.vo = sim_out.vo(1,end-n_points:end);
sim_out.vref = sim_out.vref(1,end-n_points:end);
sim_out.rmse = sum((sim_out.vo-sim_out.vref).^2)/length(sim_out.vo);


close all
figure(1)
plot(sim_out.time,sim_out.vo,'DisplayName','vo'); hold on;
plot(sim_out.time,sim_out.vref,'--','LineWidth',2,'DisplayName','vref');
legend
grid on;

%% Bode plot check
% Let  state variable x = [iL; vo]
syms Ds s

A1 = [-Rdson/L -1/L;1/Cdc -1/Ro/Cdc];
A2 = [-Rdson/L -1/L;1/Cdc -1/Ro/Cdc];
B1 = [1/L ; 0];
B2 = [0 ; 0];
U = Vg;
C = [0 , 1]; % y = vo
I2 = [1 0;0 1];
A = Ds*A1+ (1-Ds)*A2;
B = [(B1-B2)*U,Ds*B1+(1-Ds)*B2];


Gyu_s = C*inv(s*I2-A)*B
Gvod_s = Gyu_s(1,1)

D = 0.5;

Gvod = subs(Gvod_s,Ds, D)

[n_Gvod_s,d_Gvod_s] = numden(Gvod) 
n_Gvod1=double(coeffs(n_Gvod_s,'ALL'))
d_Gvod1=double((coeffs(d_Gvod_s,'ALL')))
H_Gvod1 = tf([n_Gvod1],[d_Gvod1])

figure(11)
bode(H_Gvod1)
grid on;

% Let  state variable x = [iL; vo]
%%


%% Kp KI tuning with DE
DE_out = DE(x_input,NP,CR,F,ite,@MO_DEMO_BUCK);

%%
Kp = DE_out.best_sol(1,1);
Ki = DE_out.best_sol(1,2);

temp_sim_out2 = sim('buck_demo.slx');
sim_out2.time = temp_sim_out2.sim_out.Time';
sim_out2.vo = temp_sim_out2.sim_out.Data(:,1)';
sim_out2.vref = temp_sim_out2.sim_out.Data(:,2)';

sim_out2.time = sim_out2.time(1,end-n_points:end);
sim_out2.vo = sim_out2.vo(1,end-n_points:end);
sim_out2.vref = sim_out2.vref(1,end-n_points:end);
sim_out2.rmse = sum((sim_out2.vo-sim_out2.vref).^2)/length(sim_out2.vo);

close all
figure(2)
plot(sim_out.time,sim_out.vo,'DisplayName','vo'); hold on;
plot(sim_out.time,sim_out.vref,'--','LineWidth',2,'DisplayName','vref');
plot(sim_out2.time,sim_out2.vo,'DisplayName','vo2');
legend
grid on;
%%


    
