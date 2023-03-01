%% LAB 10 -> Recursive ARX identification
% Valcauan Adina - Diana 30331/2
clc; clear; close all;

%% loaded the data
load('lab10_6.mat');

u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

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
na = 3 * n;
nb = 3 * n;

% initialize the P^(-1) matrix
P_pow_minus1 = 100 * eye(na+nb);

theta_h = recursive_arx(id, na, nb, P_pow_minus1); % using the created function

%% the final plot
sampling_T = id.Ts;

theta1 = transpose(theta_h(1:na));
theta2 = transpose(theta_h(na+1:na+nb));
A_model = [1 theta1];
B_model = [0 theta2];
C_model = [];
D_model = [];
F_model = [];

model = idpoly(A_model, B_model, C_model, D_model, F_model, 0, sampling_T);
 
figure
compare(val, model);
title('Simulated Response Comparison for RARX')

%% the function for the recursive ARX algorithm

function [theta_h] = recursive_arx(id, na, nb, P_pow_minus1)
% the id parameter is an iddata, so we need to take the arrays for input
% and output out
u_id = id.InputData;
y_id = id.OutputData;

% some initializations
n_id = length(id.u);
p = zeros(n_id, na+nb);
theta_h = zeros(na + nb, 1);

% here the phi matrix is created
% the algorithm is similar with the one from the last lab
for g_it = 1: n_id
    l = 1;
    while (l <= na) 
        if (g_it > l)
            p(g_it,l) = -y_id(g_it - l);
        end
    l = l + 1;
    end

    l = na + 1;
    while (l <= 2*na)
        if (g_it > l)
            p(g_it, l) = u_id(g_it - (l-na));
        end
        l = l + 1;
    end
end

% applying the formulas from the course/lab pdf
for g_it = 1: n_id
p_used = p(g_it,:); % I take a line at a time from the phi

error = y_id(g_it) - p_used * theta_h;

% creating 2 more auxiliary variables made writing the equation easier
% also I was able to verify the multiplications easier
ans1 = P_pow_minus1 * (p_used') * p_used * P_pow_minus1;
ans2 = 1 + p_used * P_pow_minus1 * p_used';
P_pow_minus1 = P_pow_minus1 - (ans1 / ans2);

% W = zeros(na + nb, 1);
W = P_pow_minus1 * p_used';

% this is out final theta, that I need for the graph
theta_h = theta_h + W*error; % the output of the function
end

end
