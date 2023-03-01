%% LAB 8 -> Identifying output error models using the Gauss-Newton method
% Valcauan Adina - Diana 30331/2
clc; clear; close all;

% loaded the data
load('lab8_2.mat');

u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

sampling_T = val.Ts; % sampling time

% plots for the identification and validation data
figure
subplot(211); plot(u_id, 'magenta')
title("Input identification data", "Color", 'blue')
xlabel('time'), ylabel('u');
legend('u_i_d');
subplot(212); plot(y_id, 'magenta')
title("Output identification data", "Color", 'blue')
xlabel('time'), ylabel('y');
legend('y_i_d');

figure
subplot(211); plot(u_val, 'magenta')
title("Input validation data", "Color", 'red')
xlabel('time'), ylabel('u');
legend('u_v_a_l');
subplot(212); plot(y_val, 'magenta')
title("Output validation data", "Color", 'red')
xlabel('time'), ylabel('y');
legend('y_v_a_l');

%%
% tuned for a better performance
alpha = 0.07; % stepsize -> between [0.01; 0.5]
l_iteration = 0; % iteration
lmax = 185; % [100; 250]
threshold = 1e-2; % [1e-4; 1e-1]

Nu_id = length(u_id); % the length of the identification data

f = 1; b = 1; % the initial length is 1, but it will be changed during the while loop
% initialize both the theta matrices with 1 for the beginning
theta_l = [f; b]; 
theta_lplus1 = [f; b];

% a while loop which will stop when one of the conditions is false
% for the norm, I used the "norm" function for simplicity
while (norm(theta_lplus1 - theta_l) < threshold && l_iteration <= lmax)
    % these are changed for every iteration l_iteration
    e(1) = y_id(1);
    de_f(1) = 0;
    de_b(1) = 0;
    dv = [0; 0];
    sum_dv = [0; 0]; sum_H = [0; 0];
    %dv_theta = zeros(2, N);

    for k = 2:Nu_id
    e(k) = -f*e(k-1)+y_id(k)+f*y_id(k-1)-b*u_id(k-1); % prediction error signal
    
    % derivatives
    de_f(k) = -e(k-1)-f*de_f(k-1)+y_id(k-1); 
    de_b(k) = -f*de_b(k-1)-u_id(k-1);
    end

    de_theta = [de_f; de_b]; 
    
    for k = 1:Nu_id
    sum_dv = sum_dv + e(k)*de_theta(:,k); % the sum needed for dv_theta
    sum_H = sum_H + de_theta(:,k) * transpose(de_theta(:,k)); % the sum needed for H
    end
    
    dv_theta = 2/Nu_id*sum_dv; % minimizing criterion
    H = 2/Nu_id*sum_H;
    
    
    theta_lplus1 = theta_l - alpha * inv(H) * dv_theta;
    theta_l = theta_lplus1;
    f = theta_l(1); % f gets the first element from the theta matrix
    b = theta_l(2); % b gets the last element from the theta matrix

    l_iteration = l_iteration+1;
end

%% final plot
figure
% we use idpoly and compare
% A,C,D,E are initialized with 1
model_armax = idpoly(1, [0 b], 1, 1, [1 f], 0, sampling_T);
compare(model_armax, val);