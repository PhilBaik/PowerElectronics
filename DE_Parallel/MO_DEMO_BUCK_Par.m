function f = MO_DEMO_BUCK_Par(x)
    global Fs Fsamp Ts Tsamp Vg Rdson Tend 

    mdl = 'buck_converter_DE_Parallel';
    isModelOpen = bdIsLoaded(mdl);
    open_system(mdl);
    
    numSims = length(x(:,1));

for i = 1:1:numSims
    in_temp(i) = Simulink.SimulationInput(mdl);
    in(i) = setVariable(in_temp(i),'Kp',x(i,1));
    in(i) = setVariable(in(i),'Ki',x(i,2));
    in(i) = setBlockParameter(in(i),[mdl '/Buck_converter'], 'L', num2str(x(i,3)));
end

    out = parsim(in, 'ShowProgress','off','TransferBaseWorkspaceVariables','on');
    

for i = 1:1:numSims
        simOut = out(i);
        vo = simOut.logsout.get('vo').Values;
        vref = simOut.logsout.get('vref').Values;
        duty = simOut.logsout.get('Duty').Values;
        f(i,1)=sum((vo.Data(:,1)-vref.Data(:,1)).^2)/length(vo.Data(1,:))+1e8*x(i,3);
end
    % plot(vo);hold on;plot(vref)
    % length(vo.Data(:,1))-length(vref.Data(:,1))
    % disp('hi')
end