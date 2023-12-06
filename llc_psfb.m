%% PS FB
C = 1e-11;
C_r = 0.001;
A = 10;

%% LLC
clc
clear

Ts_Control = 2.5e-5 %% sample time
Cs_fet = inf        %% snubber capacitance of ideal switch 
Rs_fet = 1e8  %% snubber resistance of ideal switch
Vf = 0.6            %% Body diode drop 

Vref = 24           %% 24

Kp = 100
Ki = 15e6
Td = 1e-9
