%% LAB 9 -> Instrumental variable methods
% Valcauan Adina - Diana 30331/2
clc; clear; close all;

% loaded the data
load('lab9_5.mat');

n_optimal = n; % given in the .mat file as the optimal value
na_optim = n_optimal; nb_optim = n_optimal; % they take the n value
nk_delay = 1; % the delay

u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

%% plots for the identification and validation data
% identification
figure
subplot(211); plot(u_id, 'magenta')
title("Input identification data", "Color", 'blue')
xlabel('time'), ylabel('u');
legend('u_i_d');
subplot(212); plot(y_id, 'magenta')
title("Output identification data", "Color", 'blue')
xlabel('time'), ylabel('y');
legend('y_i_d');

% validation
figure
subplot(211); plot(u_val, 'magenta')
title("Input validation data", "Color", 'red')
xlabel('time'), ylabel('u');
legend('u_v_a_l');
subplot(212); plot(y_val, 'magenta')
title("Output validation data", "Color", 'red')
xlabel('time'), ylabel('y');
legend('y_v_a_l');

%% the part with arx
figure
model_arx = arx(id, [na_optim nb_optim nk_delay]);
compare(model_arx, val);
title('Simulated response comparison only with ARX');

y_h = sim(model_arx, u_id); % to optain the y hat needed for the z matrix

%% z(k)
n_h = length(y_h);

z_matrix = zeros(n_h, 2*na_optim); % preallocate the matrix

% the algorithm is similar with the one in lab 6, but with the some
% modification that were needed

% first part with the output given by y hat (used sim function above)
for g_it = 1:n_h-1
    l = 1;
    while (l <= na_optim)
        if (g_it > l) 
            z_matrix(g_it, l) = -y_h(g_it - l);
        end
        l = l + 1;
    end
end

% second part with the identification input
for g_it = 1:n_h-1
    l = na_optim + 1;
    while (l <= 2*na_optim)
        if (g_it > l) 
            z_matrix(g_it, l) = u_id(g_it - (l-na_optim));
        end
        l = l + 1;
    end
end

%% phi
n_id = length(u_id);

p = zeros(n_h, 2*na_optim); % preallocate the matrix

% the algorithm is similar with the one in lab 6, but with the some
% modification that were needed

% first part with the output given by the output for identification
for g_it = 1:n_id-1
    l = 1;
    while (l <= na_optim) 
        if (g_it > l)
            p(g_it,l) = -y_id(g_it - l);
        end
    l = l + 1;
    end
end

% second part with the identification input
for g_it = 1:n_id-1
    l = na_optim + 1;
    while (l <= 2*na_optim)
        if (g_it > l)
            p(g_it, l) = u_id(g_it - (l-na_optim));
        end
        l = l + 1;
    end
end

%% phi_tilda
% p = transpose(p);
z_matrix = transpose(z_matrix); % using the transpose to be able to do the mathematical operations correctly
p_tilda = zeros(na_optim + nb_optim, na_optim + nb_optim); % preallocate the matrix

for k = 1:n_id
    % the rows/columns used for multiplication
    z_needed = z_matrix(:,k);
    p_needed = p(k,:);

    p_tilda = p_tilda + z_needed * p_needed; % formula given in the lecture
  
end
 p_tilda = (1/n_id) * p_tilda;

%% y_tilda
% z = transpose(z);
% y_id = transpose(y_id);

y_tilda = zeros(na_optim + nb_optim, 1); % preallocate the matrix
for k = 1:n_id
    % the rows/columns used for multiplication
    z_needed = z_matrix(:,k);
    y_needed = y_id(k);

    y_tilda = y_tilda + z_needed * y_needed; % formula given in the lecture
end
y_tilda = (1/n_id) * y_tilda;


%% theta
theta_needed = p_tilda \ y_tilda;

%% iv model using idploy
sampling_T = id.Ts;

theta1 = transpose(theta_needed(1:na_optim)); % used for A parameter
theta2 = transpose(theta_needed(na_optim+1:2*na_optim)); % used for B parameter
A_model = [1 theta1];
B_model = [0 theta2];
C_model = []; 
D_model = [];
F_model = [];

model_iv = idpoly(A_model, B_model, C_model, D_model, F_model, 0, sampling_T); % obtaining the model

%% the final plot
figure
% graphs with the 2 models compared with the validation data
compare(val, model_iv, model_arx); 
title('Final simulated response comparison with both ARX and IV');
