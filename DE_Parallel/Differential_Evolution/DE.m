function out = DE(x,NP,CR,F,ite,MO)
    
    out.N = length(x.max);
    for i = 1:1:NP
            out.population(i,:,1) = rand(1,out.N).*(x.max-x.min)+x.min;
    end
    out.y(:,:,1) = MO(out.population(:,:,1));

    for k = 1:1:ite
        for i = 1:1:NP
            temp1 = randi([1 NP]);
            temp2 = randi([1 NP]);
            temp3 = randi([1 NP]);
        
            out.donor(i,:) = out.population(temp1,:,k) + F*(out.population(temp2,:,k)-out.population(temp3,:,k));

            out.CR_comparison_temp = rand([1, out.N])<CR;
            out.trial(i,:) = out.donor(i,:).*out.CR_comparison_temp +out.population(i,:,k).*~out.CR_comparison_temp;

        end
        out.trial;

        % By using Matrix
        MO_trial = MO(out.trial(:,:));
        MO_population = MO(out.population(:,:,k));

        out.population(:,:,k+1) = out.trial.*(MO_trial<MO_population) + out.population(:,:,k).*(MO_trial>=MO_population);
        out.y(:,:,k+1) = MO_trial.*(MO_trial<MO_population) + MO_population.*(MO_trial>=MO_population);

        % By using For loop
        % for m = 1:1:NP
        %     if MO_trial(m,1)<MO_population(m,1)
        %         out.population(m,:,k+1) = out.trial(m,:);
        %         out.y(m,:,k+1) = MO_trial(m,1);
        %     else
        %         out.population(m,:,k+1) = out.population(m,:,k);
        %         out.y(m,:,k+1) = MO_population(m,1);
        %     end
        % end
        [value out.best_index] = min(out.y(:,1,end));
        out.best_sol = out.population(out.best_index,:,end);


        disp_name = sprintf('iteration: %d',k);
        disp(disp_name);
    end

    

end