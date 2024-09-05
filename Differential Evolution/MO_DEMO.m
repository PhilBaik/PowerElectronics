function [f partial_f] = MO_DEMO(varargin)
    if length(varargin)==1
        x = varargin{1};
        x1 = x(:,1);
        x2 = x(:,2);
    else
        x1 = varargin{1};
        x2 = varargin{2};
    end
    
    f = (x1 - 2).^2 + (x2 - 3).^2 + sin(5 * x1) .* sin(5 * x2);
    partial_f = f;
end

%% DE
% function out = DE(x,NP,CR,F,ite,MO)
% 
%     out.N = length(x.max);
%     for i = 1:1:NP
%             out.population(i,:,1) = rand(1,out.N).*(x.max-x.min)+x.min;
%             out.y(i,:,1) = MO(out.population(i,:,1));
%     end
% 
%     for k = 1:1:ite
%         for i = 1:1:NP
%             temp1 = randi([1 NP]);
%             temp2 = randi([1 NP]);
%             temp3 = randi([1 NP]);
% 
%             out.donor(i,:) = out.population(temp1,:,k) + F*(out.population(temp2,:,k)-out.population(temp3,:,k));
% 
%             out.CR_comparison_temp = rand([1, out.N])<CR;
%             out.trial(i,:) = out.donor(i,:).*out.CR_comparison_temp +out.population(i,:,k).*~out.CR_comparison_temp;
%             MO_trial = MO(out.trial(:,:));
%             MO_population = MO(out.population(:,:,k));
% 
%             for m = 1:1:NP
%                 if MO_trial(m,1)<MO_population(m,1)
%                     out.population(i,:,k+1) = out.trial(i,:);
%                     out.y(i,:,k+1) = MO_trial;
%                 else
%                     out.population(i,:,k+1) = out.population(i,:,k);
%                     out.y(i,:,k+1) = MO_population;
%                 end
%             end
%             disp_name = sprintf('i: %d',i);
%             disp(disp_name);
%         end
%         disp_name = sprintf('iteration: %d',k);
%         disp(disp_name);
%     end
% 
%     [value out.best_index] = min(out.y(:,1,end));
%     out.best_sol = out.population(out.best_index,:,end);
% 
% end