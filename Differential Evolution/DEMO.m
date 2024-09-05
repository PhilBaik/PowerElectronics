%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Title: Differential Evolution (DEMO)
%%%%%%%%%%%% Writer: Hyeongmeen Baik
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% Initialization
xi1_min = -2
xi1_max = 5
xi2_min = -2
xi2_max = 5

NP = 10     %Population size: 10N (N - number of dimension)

x_input.min = [xi1_min,xi2_min];
x_input.max = [xi1_max,xi2_max];
xi1 = linspace(xi1_min,xi1_max,300);
xi2 = linspace(xi2_min,xi2_max,300);


CR = 0.9   %Crossover probability
F = 0.8   %Differential Weight
ite = 500;

% Meshgrid for comparison
[X1, X2] = meshgrid(xi1,xi2);
Y = MO(X1,X2);
[min_Y_col,min_col_arr] = min(Y,[],1);
[min_Y,min_col] = min(min_Y_col);
min_row = min_col_arr(1,min_col);

MO_DEMO(X1(min_row,min_col),X2(min_row,min_col))
target = [X1(min_row,min_col),X2(min_row,min_col)]
%% Differential Evolution
% Differential evolution function needs parameter input vector, NP, CR, F,
% and the number of iteration, and function handle for multi-objective
% function

DE_out = DE(x_input,NP,CR,F,ite,@MO_DEMO);
best_sol = DE_out.best_sol;
best_Y = min(min(DE_out.y));
%% Figures
close all
figure(1)
mesh(X1,X2,Y,'DisplayName','Plane');hold on
stem3(X1(min_row,min_col),X2(min_row,min_col),min_Y,'diamondr','MarkerSize',10,'LineStyle','none','DisplayName','Global minimum')
stem3(DE_out.population(:,1,1),DE_out.population(:,2,1),DE_out.y(:,1,1),'filled','LineStyle',"none",'DisplayName','initial');
stem3(DE_out.population(:,1,ite+1),DE_out.population(:,2,ite+1),DE_out.y(:,1,ite+1),'filled','MarkerSize',10,'LineStyle',"none",'DisplayName','final');
legend



    
