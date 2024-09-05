function out = DE(x,NP,CR,F,ite,MO)
    
    out.N = length(x.max);
    for i = 1:1:NP
            out.population(i,:,1) = (rand(1,out.N).*(x.max-x.min)+x.min);
    end
    out.xmin = x.min;
    out.xmax = x.max;
    [out.ytotal(:,:,1), out.ypartial(:,:,1)] = MO(out.population(:,:,1));
     

    [value out.best_index(1,1)] = min(out.ytotal(:,1,1));
    out.best_sol(:,:,1) = out.population(out.best_index(1,1),:,1);

    for k = 1:1:ite
        for i = 1:1:NP
            temp1 = randi([1 NP]);
            temp2 = randi([1 NP]);
            temp3 = randi([1 NP]);
        

            out.donor(i,:) = out.population(temp1,:,k) + F*(out.population(temp2,:,k)-out.population(temp3,:,k));
            sanity_test = ((out.donor(i,:)-out.xmin).*(out.xmax-out.donor(i,:))<0);
            sum_sanity_test = sum(sanity_test);
            while sum_sanity_test > 0
                temp1 = randi([1 NP]);
                temp2 = randi([1 NP]);
                temp3 = randi([1 NP]);
                out.donor(i,:) = out.population(temp1,:,k) + F*(out.population(temp2,:,k)-out.population(temp3,:,k));
                sanity_test = ((out.donor(i,:)-out.xmin).*(out.xmax-out.donor(i,:))<0);
                sum_sanity_test = sum(sanity_test);
            end

            out.CR_comparison_temp = rand([1, out.N])<CR;
            out.trial(i,:) = out.donor(i,:).*out.CR_comparison_temp +out.population(i,:,k).*~out.CR_comparison_temp;

        end
        out.trial;

        % By using Matrix
        [MO_trial, MO_trial_partial] = MO(out.trial(:,:));
        MO_population = out.ytotal(:,:,k);
        MO_population_partial = out.ypartial(:,:,k);

        out.population(:,:,k+1) = out.trial.*(MO_trial<MO_population) + out.population(:,:,k).*(MO_trial>=MO_population);
        out.ytotal(:,:,k+1) = MO_trial.*(MO_trial<MO_population) + MO_population.*(MO_trial>=MO_population);
        out.ypartial(:,:,k+1) = MO_trial_partial.*(MO_trial<MO_population) + MO_population_partial.*(MO_trial>=MO_population);

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
        [value out.best_index(1,k+1)] = min(out.ytotal(:,1,k+1));
        out.best_sol(:,:,k+1) = out.population(out.best_index(1,k+1),:,k+1);

        disp_name = sprintf('iteration: %d',k);
        disp(disp_name);
    end

    

end