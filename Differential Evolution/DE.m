function out = DE(x,NP,CR,F,ite,MO)
    
    out.N = length(x.max);
    for i = 1:1:NP
            out.population(i,:,1) = rand(1,out.N).*(x.max-x.min)+x.min;
            out.y(i,:,1) = MO(out.population(i,:,1));
    end

    for k = 1:1:ite
        for i = 1:1:NP
            temp1 = randi([1 NP]);
            temp2 = randi([1 NP]);
            temp3 = randi([1 NP]);
        
            out.donor(i,:) = out.population(temp1,:,k) + F*(out.population(temp2,:,k)-out.population(temp3,:,k));

            out.CR_comparison_temp = rand([1, out.N])<CR;
            out.trial(i,:) = out.donor(i,:).*out.CR_comparison_temp +out.population(i,:,k).*~out.CR_comparison_temp;
            MO_trial = MO(out.trial(i,:));
            MO_population = MO(out.population(i,:,k));

            if MO_trial<MO_population
                out.population(i,:,k+1) = out.trial(i,:);
                out.y(i,:,k+1) = MO_trial;
            else
                out.population(i,:,k+1) = out.population(i,:,k);
                out.y(i,:,k+1) = MO_population;
            end
            disp_name = sprintf('i: %d',i);
            disp(disp_name);
        end
        disp_name = sprintf('iteration: %d',k);
        disp(disp_name);
    end

    [value out.best_index] = min(out.y(:,1,end));
    out.best_sol = out.population(out.best_index,:,end);

end