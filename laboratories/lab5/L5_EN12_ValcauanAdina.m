%% Lab 5 -> Correlation Analysis
% Valcauan Adina-Diana 30331

clc; clear;

% loaded the data
load('lab5_4.mat');
u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

figure
plot(tid, u_id, 'Color', 'r'); hold on
plot(tid, y_id, 'Color', 'b'); hold off
title('Plot of the identification data'), legend('u_i_d', 'y_i_d');
xlabel('t'), ylabel('u, y');

figure
subplot(211)
plot(tval, u_val, 'Color', 'r'); hold on
title('Plot of the validation data'), legend('u_v_a_l');
xlabel('t'), ylabel('u');
subplot(212)
plot(tval, y_val, 'Color', 'b'); hold off
legend('y_v_a_l');
xlabel('t'), ylabel('u');

% using detrend if needed
u_id_d = detrend(u_id);
y_id_d = detrend(y_id);

figure
plot(tid, u_id_d, 'Color', 'r'); hold on
plot(tid, y_id_d, 'Color', 'b'); hold off
title('Plot of the identification after detrend'), legend('u_i_d', 'y_i_d');
xlabel('t'), ylabel('u, y');

% r_u and r_uy -> covariance functions
N = length(y_id);

for tau = 1 : N
    % initialize the sum and the array element with 0
    s1 = 0; s2 = 0;
    r_u_array (tau) = 0; r_uy_array(tau) = 0;

    for k = 0 : (N - tau)
    s1 = s1 + (y_id(k + tau) * u_id_d(k + 1));
    s2 = s2 + (u_id_d(k + tau) * u_id_d(k + 1));
    
    r_u_array(tau) = r_u_array(tau) + (1/N) * s2;
    r_uy_array(tau) = r_uy_array(tau) + (1/N) * s1;
    end
end


% R
T = length(u_id);

figure
stem(imp);
M = 66;

for i = 1:T
    for j = 1:M
        a_matrix(i,j) = r_u_array(abs(j-i)+1);
    end
end

% FIR model
h_h = a_matrix\transpose(r_uy_array);

% y hat for identification and validation
y_h_id = conv(h_h, u_id_d);
y_h_val = conv(h_h, u_val);

% MSE for identification and validation
N_id = length(y_id_d);
y_h_id_cut = y_h_id(1:2500); % the y hat has too many values, so we cut it
mse_id = (1 / N_id) * sum((y_h_id_cut - y_id_d).^2);

N_val = length(y_val);
y_h_val_cut = y_h_val(1:250); 
mse_val = (1 / N_val) * sum((y_h_val_cut - y_val).^2);

% plots
figure
plot(tid, y_id_d); hold on
plot(tid, y_h_id_cut);
title('Final plot: identification');
legend('System', 'FIR Model');


figure
plot(tval, y_val); hold on
plot(tval, y_h_val_cut);
title('Final plot: validation');
legend('System', 'FIR Model');
