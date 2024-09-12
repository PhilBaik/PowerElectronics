restoredefaultpath
delete(gcp('nocreate'));
clc;
clear all;
close all;

%% Case 1
nummodel = 1;
numcontrol = 32;
numinputsim = 16;

numSims = nummodel*numcontrol*numinputsim;

mdl = 'buck_converter_data_gen';
isModelOpen = bdIsLoaded(mdl);
open_system(mdl);

%min max parameters
Fs_min = 1e3;
Fs_max = 1e6;
Rdson_min = 1e-4;
Rdson_max = 1e-2;
Rdon_min = 1e-4;
Rdon_max = 1e-2;
Kp_min = 0;
Kp_max = 1;
Vref_min = 10;
Vref_max = 90;
Ro_min = 1;
Ro_max = 100;
L_min = 0;
L_max = 1e-3;
Co_min = 0;
Co_max = 1e-8;
Vg_min = 10;
Vg_max = 3300;
Vf_min = 0;
Vf_max = 2;

%% Test
i=1;
Fs(1,i) = rand_gen(Fs_min,Fs_max);
Fsamp(1,i) = Fs(1,i)*100;
Ts(1,i) = 1/Fs(1,i);
Tsamp(1,i) = 1/Fsamp(1,i);

Rdson(1,i) = rand_gen(Rdson_min,Rdson_max);
Rdon(1,i) = rand_gen(Rdon_min,Rdon_max);

Kp(1,i) = rand_gen(Kp_min,Kp_max);
Ki(1,i) = rand_gen(Kp_min,Kp_max);

Tend(1,i) = Ts(1,i)*1000;

Ro(1,i) = rand_gen(Ro_min,Ro_max);
Ttrig1(1,i) = Tend(1,i)/3;
Ttrig2(1,i) = Tend(1,i)*2/3;

Vg(1,i) = rand_gen(Vg_min,Vg_max);
Vf(1,i) = rand_gen(Vf_min,Vf_max);

Vref1(1,i) = rand_gen(Vref_min,Vref_max);
Vref2(1,i) = rand_gen(Vref_min,Vref_max);
Vref3(1,i) = rand_gen(Vref_min,Vref_max);
L(1,i) = rand_gen(L_min,L_max);
Co(1,i) = rand_gen(Co_min,Co_max);
Tend(1,i) = Ts(1,i)*1000;
%%

%%
for l = 1:1:nummodel
    Fs(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol)= rand_gen(Fs_min,Fs_max);
    Fsamp(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = Fs(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol)*500;
    Ts(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = 1/Fs(1,(l-1)*numinputsim*numcontrol+1);
    Tsamp(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = 1/Fsamp(1,(l-1)*numinputsim*numcontrol+1);
    Tend(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = Ts(1,(l-1)*numinputsim*numcontrol+1)*1000;
    Ttrig1(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = Tend(1,(l-1)*numinputsim*numcontrol+1)/3;
    Ttrig2(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = Tend(1,(l-1)*numinputsim*numcontrol+1)*2/3;

    Rdon(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = rand_gen(Rdon_min,Rdon_max);
    Rdson(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = rand_gen(Rdson_min,Rdson_max);
    Vg(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = rand_gen(Vg_min,Vg_max);
    Vf(1,(l-1)*numinputsim*numcontrol+1:l*numinputsim*numcontrol) = rand_gen(Vf_min,Vf_max);
    for m = 1:1:numcontrol
        Kp(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+1:(l-1)*numinputsim*numcontrol+m*numinputsim) = rand_gen(Kp_min,Kp_max);
        Ki(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+1:(l-1)*numinputsim*numcontrol+m*numinputsim) = rand_gen(Kp_min,Kp_max);

        for n = 1:1:numinputsim
  

        Ro(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+n) = rand_gen(Ro_min,Ro_max);
        
        Vref1(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+n) = rand_gen(Vref_min,Vref_max);
        Vref2(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+n) = rand_gen(Vref_min,Vref_max);
        Vref3(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+n) = rand_gen(Vref_min,Vref_max);
        L(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+n) = rand_gen(L_min,L_max);
        Co(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+n) = rand_gen(Co_min,Co_max);
        end
    end
end

%%
for i = 1:1:numSims
    in_temp(i) = Simulink.SimulationInput(mdl);
    in_temp(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'R', num2str(Ro(1,i)));
    in_temp(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'Rdon', num2str(Rdon(1,i)));
    in_temp(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'Rdson', num2str(Rdson(1,i)));
    in_temp(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'Vg', num2str(Vg(1,i)));
    in_temp(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'Vf', num2str(Vf(1,i)));
    in_temp(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'C', num2str(Co(1,i)));
    in_temp(i) = setBlockParameter(in_temp(i),[mdl '/Buck_converter'], 'L', num2str(L(1,i)));
    in_temp(i) = setVariable(in_temp(i),'Vref1', Vref1(1,i));
    in_temp(i) = setVariable(in_temp(i),'Vref2', Vref2(1,i));
    in_temp(i) = setVariable(in_temp(i),'Vref3', Vref3(1,i));
    in_temp(i) = setVariable(in_temp(i),'Ttrig1', Ttrig1(1,i));
    in_temp(i) = setVariable(in_temp(i),'Ttrig2', Ttrig2(1,i));
    in_temp(i) = setVariable(in_temp(i),'Kp', Kp(1,i));
    in_temp(i) = setVariable(in_temp(i),'Ki', Ki(1,i));
    in_temp(i) = setVariable(in_temp(i),'Ts', Ts(1,i));
    in_temp(i) = setVariable(in_temp(i),'Tsamp', Tsamp(1,i));
    in_temp(i) = setVariable(in_temp(i),'Tend', Tend(1,i));
    in(i) = in_temp(i);
end
%%
out_sim = parsim(in, 'ShowProgress','on','TransferBaseWorkspaceVariables','off');
delete(gcp('nocreate'));
if(~isModelOpen)
    close_system(mdl, 0);
end

%%
simData = cell(1,numSims);
%%
for i = 1:1:numSims
    simData{i}.Variables = in(i).Variables;
    simData{i}.BlockParameters = in(i).BlockParameters;
    simData{i}.time = out_sim(i).logsout.get('vo').Values.Time';
    simData{i}.vo = out_sim(i).logsout.get('vo').Values.Data(:,1)';
    simData{i}.vref = out_sim(i).logsout.get('vref').Values.Data(:,1)';
    simData{i}.gate = out_sim(i).logsout.get('Gate').Values.Data(:,1)';
    simData{i}.io = out_sim(i).logsout.get('io').Values.Data(:,1)';
end
%%
i=1;
foldername = sprintf('%s',date());
mkdir(foldername);
filename1 = sprintf('%s/ver%d_1_1.mat',date(),i);
while exist(filename1)
    i = i+1;
    filename1 = sprintf('%s/ver%d_1_1.mat',date(),i);
end

for l = 1:1:nummodel
    for m = 1:1:numcontrol
        filename_temp = sprintf('%s/ver%d_%d_%d.mat',date(),i,l,m);
        temp_Data = simData(1,(l-1)*numinputsim*numcontrol+(m-1)*numinputsim+1:(l-1)*numinputsim*numcontrol+m*numinputsim);
        save(filename_temp,'temp_Data')
    end
end


%%
close all

figure(1)
i=1;
simOut = temp_Data{1,1};
% 'ts' is a MATLAB 'timeseries' object that stores the time and
% data values for the logged 'vertical_disp' signal.
% Use the 'plot' method of the object to plot the data against the
% time.
plot(simOut.time,simOut.vo,'DisplayName','vo, i=1');
hold on
plot(simOut.time,simOut.vref,'DisplayName','vref, i=1');
title('i = 1')

figure(2)
simOut = temp_Data{1,end};
% 'ts' is a MATLAB 'timeseries' object that stores the time and
% data values for the logged 'vertical_disp' signal.
% Use the 'plot' method of the object to plot the data against the
% time.
plot(simOut.time,simOut.vo,'DisplayName','vo, i=1');
hold on
plot(simOut.time,simOut.vref,'DisplayName','vref, i=1');
title('i = end')
%% local func
function y= rand_gen(min_x, max_x)
    y=rand(1,1).*(max_x-min_x)+min_x;
end