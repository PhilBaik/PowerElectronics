clc;
clear;

for ite = 1:1:20
    run('data_generator_case1_1.m')
    run('data_generator_case1_2.m')
end


% hi1 = load('ver2_1_1.mat','temp_Data')
% hi2 = load('ver2_1_2.mat','temp_Data')