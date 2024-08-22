function f = MO_DEMO_BUCK(x)
    global Fs Fsamp Ts Tsamp Vg Rdson Tend Kp Ki
    Kp = x(1,1);
    Ki = x(1,2);

    temp_sim_out = sim('buck_demo.slx');
    sim_out.time = temp_sim_out.sim_out.Time';
    sim_out.vo = temp_sim_out.sim_out.Data(:,1)';
    sim_out.vref = temp_sim_out.sim_out.Data(:,2)';
    
    n_cycles = 100;
    n_points = n_cycles*Ts/Tsamp;
    
    sim_out.time = sim_out.time(1,end-n_points:end);
    sim_out.vo = sim_out.vo(1,end-n_points:end);
    sim_out.vref = sim_out.vref(1,end-n_points:end);

    f = sum((sim_out.vo-sim_out.vref).^2)/length(sim_out.vo);

end