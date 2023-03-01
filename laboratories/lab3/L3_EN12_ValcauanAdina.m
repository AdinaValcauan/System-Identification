%% Lab 3 -> Transient Analysis of Step Responses
% Valcauan Adina-Diana 30331

%% first order system
clc; clear;

% loaded the data
load('lab3_order1_5.mat');
u_data_1st = data.InputData;
y_data_1st = data.OutputData;

% first 100 values
u_id_1st = u_data_1st(1:100,1);
y_id_1st = y_data_1st(1:100,1);

length_time_1st = length(t);
floor_time_length_1st = floor(length_time_1st/5);
t_id_1st = t(1:floor_time_length_1st);


plot(t_id_1st, y_id_1st, 'LineWidth', 2, 'Color', 'r'); hold on
plot(t_id_1st, u_id_1st, 'LineWidth', 2, 'Color', 'b');

% reading from the graphic inital/steady-state values for i/o
u0 = 0; u_ss = 4;
y0 = 0.02; y_ss = 1.98;

% calculating k (gain)
k = (y_ss - y0) / (u_ss - u0)

y_calc = y0 + 0.6322*(y_ss-y0)
T = 3.43; 

% building the first order transfer function
H1 = tf(k, [T 1])

% MSE
y_lsim = lsim(H1,u_data_1st, t);
mse = (1 / length_time_1st) .* (y_data_1st - y_lsim).^2


% final plot
figure
plot(t, y_lsim); hold on
plot(t, y_data_1st);

%% second order system
clc; clear;

% loaded the data
load('lab3_order2_5.mat');
u_data_2nd = data.InputData;
y_data_2nd = data.OutputData;

% first 100 values
u_id_2nd = u_data_2nd(1:100,1);
y_id_2nd = y_data_2nd(1:100,1);

length_time_2nd = length(t);
floor_time_length_2nd = floor(length_time_2nd/5);
t_id_2nd = t(1:floor_time_length_2nd);


plot(t_id_2nd, y_id_2nd, 'LineWidth', 2, 'Color', 'r'); hold on
plot(t_id_2nd,u_id_2nd, 'LineWidth', 2, 'Color', 'b');

% reading from the graphic inital/steady-state values for i/o
u0 = 0; u_ss = 0.5; 
y0 = -0.06; y_ss = 0.97;

k = (y_ss - y0) / (u_ss - u0) % the gain

% reading the time/ratio for first peak and first valley
t1= 0.78; y_t1 = 1.58; % first peak
t2 = 1.48; y_t2 = 0.76; % first valley

M = (y_ss - y_t2) / (y_t1 - y_ss) % the overshoot
T0 = 2 * (t2 - t1) % oscillation period
tita = (log(1/M)) / (sqrt(pi*pi + log(M)*log(M))) % the damping
wn = (2*pi) / (T0 * sqrt(1-tita*tita)) % the natural frequency

% building the second order transfer function
H2 = tf(k*wn*wn, [1 2*tita*wn wn*wn])

% MSE
y_lsim = lsim(H2, u_data_2nd, t);
mse = (1 / length_time_2nd) .* (y_data_2nd - y_lsim).^2


% final plot
figure
plot(t, y_lsim); hold on
plot(t, y_data_2nd);

