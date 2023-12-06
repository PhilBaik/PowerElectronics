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

%%

%%%%% Quality factor = sqrt(Lr/Cr)/Rac
%%%%% where Rac = 8/pi^2/n^2/Ro
%%%%% m = (Lr+Lm)/Lr Ratio of total primary inductance to resonant inductance
%%%%% M = voltage gain of resonant tank or K

m = 6
Fx = 0.01:0.0001:100;


K_2 = K(0.2,m,Fx);
K_4 = K(0.3,m,Fx);
K_6 = K(0.5,m,Fx);
K_8 = K(0.7,m,Fx);
K_10 = K(1,m,Fx);
K_50 = K(5,m,Fx);

Lr = 2.2e-6
Lm = 12.2e-6
n = 12                  %%Np/Ns
Cr = 0.94e-6    
fs_min = 50e3   
fr = 110e3              %%1/(2pi sqrt(LrCr))
m = (Lr+Lm)/Lr          %%Ratio of total primary inductance to resonant inductance
M_nom = 1

%%%%% Specification
%%%%% Vo = 400
%%%%% Vin = 33 nominal
%%%%% Vin_min = 25
%%%%% Vin_max = 36
%%%%% Output power = 250
%%%%% Resonant frequency 110k
%%%%% reasonable m 6-10
%%%%% Fx = fs/fr

%%%%% Quality factor = sqrt(Lr/Cr)/Rac
%%%%% where Rac = 8/pi^2*n^2*Ro
%%%%% m = (Lr+Lm)/Lr Ratio of total primary inductance to resonant inductance
%%%%% M = voltage gain of resonant tank or part of K
Vin = 33
Vo = 400
n = Vin/Vo
Vin_min = 25
Vin_max = 36
%%%%% n*M = total gain
M_min = Vo/Vin_max*n
M_max = Vo/Vin_min*n
Po = 250

Qmax = 0.4
m = 6.3

K1 = K(Qmax,m,Fx);
diff_K1 = diff(K1);
zero_diff_K1 = diff_K1<0;

for i=1:1:length(zero_diff_K1)
    if zero_diff_K1(1,i) == 1
        index = i;
        break;
    end
end

Fx_min = Fx(1,index)
K_max = K1(1,index)

fr = 110e3
fs_min = Fx_min*fr

if(K_max > M_max)
    disp('fine enough!')
end

%%%%% full load
%%%%% where Rac = 8/pi^2*n^2*Ro
Rac_min = 8/pi^2*n^2*Vo^2/Po
Lr = Rac_min*Qmax/fr/2/pi
Cr = 1/(2*pi*Rac_min*Qmax*fr)
Lm = Lr*m-Lr
Ro = Rac_min/8*pi^2/n^2
%% Simulink


%% plot
figure(1)
semilogx(Fx,K_2,'Displayname','Q = 0.2')
title('m = 6')
grid on
hold on
semilogx(Fx,K_4,'Displayname','Q = 0.4')
semilogx(Fx,K_6,'Displayname','Q = 0.6')
semilogx(Fx,K_8,'Displayname','Q = 0.8')
semilogx(Fx,K_10,'Displayname','Q = 1')
semilogx(Fx,K_50,'Displayname','Q = 5')
legend 

figure(2)
semilogx(Fx,K1)
grid on;

%%
function out = K(Q,m,Fx)
    out = Fx.^2*(m-1)./sqrt( (m*(Fx.^2)-1).^2 + (Fx.^2).*((Fx.^2-1).^2)*((m-1)^2)*Q^2 );
end

