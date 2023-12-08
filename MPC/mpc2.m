clc
clear

I = 0.01:0.01:30

Lo = 3e-1
Ls = 0.1e-6
phi_s = 10e-1

phi = Ls*I + tanh((Lo-Ls)/phi_s*I);
Is = 1
a1 = Lo;
a2 = 3*phi_s/Is^2 - 2*Lo/Is - Ls/Is;
a3 = (Ls+Lo)/Is^2 - 2*phi_s/Is^3;

for i = I
    if i<Is
        phi2 = a1*I + a2*I.^2+a3*I.^3;
        disp('hi')
    else
        phi2 = phi_s+ Ls.*(I-Is);
        disp('no')
    end
end



I_tot = horzcat(I,-I);
phi2 = horzcat(phi2,-phi2);
I_tot = sort(I_tot);
phi2 = sort(phi2);

flux_current = horzcat(phi',I');
flux_current2 = horzcat(phi2',I_tot');
plot(I_tot,phi2) 