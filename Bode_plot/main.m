%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Title: Bode Plot
%%%%%%%%%%%% Writer: Hyeongmeen Baik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear
close all

%% Initialization
global Fs Fsamp Ts Tsamp Vg Rdson Tend Ttrig L Cdc Ro


NP =  40;    %Population size: 10N (N - number of dimension)

Kp = 0.04;
Ki = 1;

CR = 0.9;   %Crossover probability
F = 0.2;    %Differential Weight
ite = 10;

Fs = 1e5;
Fsamp = 1e6;
Ts = 1/Fs;
Tsamp = 1/Fsamp;
Rdson = 0.01;
Tend = 0.1;
Vg = 100;
Ro = 1;
Ttrig = Tend/2;

Vref1 = 40;
Vref2 = 60;

L = 1e-3;
Cdc = 1e-4;

Kp_init = 1;
Ki_init = 1;

%% Bode plot check
% Let  state variable x = [iL; vo]
syms Ds s
Rdson = 1e-3;

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

D = 0.3;

Gvod = subs(Gvod_s,Ds, D)

[n_Gvod_s,d_Gvod_s] = numden(Gvod) 
n_Gvod1=double(coeffs(n_Gvod_s,'ALL'))
d_Gvod1=double((coeffs(d_Gvod_s,'ALL')))
H_Gvod1 = tf([n_Gvod1],[d_Gvod1])
H_Gvod1_z = c2d(H_Gvod1,Tsamp,'ZOH')




ff = logspace(2,log10(Fsamp)-1,10);


%%

sim("Buck_bode.slx");
data_high = logsOut.get('Data').Values;
sys_estim_high = frestimate(data_high,ff*2*pi,'rad/s');

%%

Fs = 1e4;
Ts = 1/Fs;
sim("Buck_bode.slx");
data_low = logsOut.get('Data').Values;
sys_estim_low = frestimate(data_low,ff*2*pi,'rad/s');

%%
figure(21)
bode(H_Gvod1);hold on;
bode(H_Gvod1_z);
bode(sys_estim_high);
bode(sys_estim_low);
legend('s-theory','z-theory','fs = 1e5','fs =1e4')
grid on;

figure(22)
bode(sys_estim_high);hold on;
bode(sys_estim_low);