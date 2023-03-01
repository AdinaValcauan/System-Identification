%% Lab 4 -> Transient Analysis of Impulse Responses
% Valcauan Adina-Diana 30331

%% first order system
clc; clear;

% loaded the data
load('lab4_order1_4.mat');
u_data_1st = data.InputData;
y_data_1st = data.OutputData;

figure
plot(t, u_data_1st, 'Color', 'r'); hold on
plot(t, y_data_1st, 'Color', 'b'); hold off
title('First order: Plot of the data'), legend('u_d_a_t_a', 'y_d_a_t_a');
xlabel('t'), ylabel('u, y');

% the first impulse used for identification
u_id_1st = u_data_1st(16:77,1);
y_id_1st = y_data_1st(16:77,1);
t_id_1st = t(16:77);

figure
plot(t_id_1st, u_id_1st, 'Color', 'r'); hold on
plot(t_id_1st, y_id_1st, 'Color', 'b'); 
title('First order: Identification'), legend('u_i_d', 'y_i_d');
xlabel('t'), ylabel('u, y');

% the last 2 impulses used for validation
u_val_1st = u_data_1st(78:330,1);
y_val_1st = y_data_1st(78:330,1);
t_val_1st = t(78:330);
figure
plot(t_val_1st, u_val_1st, 'Color', 'r'); hold on
plot(t_val_1st, y_val_1st, 'Color', 'b'); hold off
title('First order: Validation'), legend('u_v_a_l', 'y_v_a_l');
xlabel('t'), ylabel('u, y');

% the inital/steady-state values for i/o
u0_1st = 0; u_ss_1st = 1.5;
y0_1st = 1.47; y_ss_1st = 1.57;

% the gain
k = y_ss_1st / u_ss_1st;

% finding the time constant
y_max = 3.42;
y_t2 = y0_1st + 0.368*(y_max - y0_1st);

t1 = 1.2; t2 = 1.8;
T = t2 - t1;

% the matrices
A = -1/T;
B = k/T;
C = 1;
D = 0;

H_ss_1st = ss(A,B,C,D)

y1 = lsim(H_ss_1st, u_val_1st, t_val_1st, y_ss_1st);

% the MSE for the validation
mse_val_1st = (1 / length(y_val_1st)) .* sum((y_val_1st - y1).^2);

% the final plot
figure
plot(t_val_1st, y1); hold on
plot(t_val_1st, y_val_1st); 
title('First oder: Final plot');
xlabel('t');


%% second order system
clc; clear;

% loaded the data
load('lab4_order2_4.mat');
u_data_2nd = data.InputData;
y_data_2nd = data.OutputData;

figure
plot(t, u_data_2nd, 'Color', 'r'); hold on
plot(t, y_data_2nd, 'Color', 'b'); hold off
title('Second order: Plot of the data'), legend('u_d_a_t_a', 'y_d_a_t_a');
xlabel('t'), ylabel('u, y');

% the first impulse used for identification
u_id_2nd = u_data_2nd(15:86,1);
y_id_2nd = y_data_2nd(15:86,1);
t_id_2nd = t(15:86);

figure
plot(t_id_2nd, u_id_2nd, 'Color', 'r'); hold on
plot(t_id_2nd, y_id_2nd, 'Color', 'b'); hold off
title('Second order: Identification'), legend('u_i_d', 'y_i_d');
xlabel('t'), ylabel('u, y');

% the last 2 impulses used for validation
u_val_2nd= u_data_2nd(86:330,1);
y_val_2nd = y_data_2nd(86:330,1);
t_val_2nd = t(86:330);

figure
plot(t_val_2nd, u_val_2nd, 'Color', 'r'); hold on
plot(t_val_2nd, y_val_2nd, 'Color', 'b'); hold off
title('Second order: Validation'), legend('u_v_a_l', 'y_v_a_l');
xlabel('t'), ylabel('u, y');

% the inital/steady-state values for i/o
u0_ = 0; u_ss_2nd = 1;
y0_2nd = 0.252; y_ss_2nd = 0.254;

k = y_ss_2nd / u_ss_2nd; % the gain

% reading the time/ratio for first peak and first valley
t1= 1.2; y_t1 = 0.93; % first peak
t2 = 1.95; y_t2 = 0.13; % first valley

% reading from the graph the needed values
t00 = 0.96; t01 = 1.69; t02 = 2.4; % the time values
y00 = 0.25; y01 = 0.255; y02 = 0.24; % the y values
k00 = 17; k01 = 40; k02 = 62; % the indexes

Ts = 2 * (t2 - t1);
A_plus = Ts * sum(y_id_2nd(k00:k01) - y0_2nd); % for the first peak
A_minus = Ts * sum(y0_2nd - y_id_2nd(k01:k02)); % for the first valey

M = A_minus / A_plus; % the overshoot
tita =  (log(1/M)) / (sqrt(pi*pi + log(M)*log(M))); % the damping

T0 = 2 * (t2 - t1); % oscillation period
wn = (2*pi) / (T0 * sqrt(1-tita*tita)); % the natural frequency

% the matrices
A = [0 1; -wn^2 -2*tita*wn];
B = [0; k*(wn^2)];
C = [1 0];
D = 0;

H_ss_2nd = ss(A,B,C,D)

y2 = lsim(H_ss_2nd, u_val_2nd, t_val_2nd); 

% the MSE for the validation
mse_val_2nd = (1 / length(y_val_2nd)) .* sum((y_val_2nd - y2).^2);

% the final plot
figure
plot(t_val_2nd, y2); hold on
plot(t_val_2nd, y_val_2nd); 
title('Second oder: Final plot');
xlabel('t');
