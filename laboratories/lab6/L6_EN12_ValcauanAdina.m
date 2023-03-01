%% LAB 6 -> ARX identification
% Valcauan Adina - Diana 30331/2

clc; clear;

% loaded the data
load('lab6_1.mat');
u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

% plot for identification
figure
subplot(2, 1, 1)
plot(u_id, 'Color', 'r'); hold on
title('Plot of the identification data'), legend('u_i_d')
xlabel('time'), ylabel('u')
subplot(2, 1, 2)
plot(y_id, 'Color', 'b'); hold off
legend('y_i_d');
xlabel('time'), ylabel('y');

% plot for validation
figure
subplot(2, 1, 1)
plot(u_val, 'Color', 'r'); hold on
title('Plot of the validation data'), legend('u_v_a_l');
xlabel('time'), ylabel('u');
subplot(2, 1, 2)
plot(y_val, 'Color', 'b'); hold off
legend('y_v_a_l'); 
xlabel('time'), ylabel('y');
 
%% finding na and nb
na_search = 1:20; nb_search = 1:20; nk_f = 1;
nn_f = struc(na_search, nb_search, nk_f);

v_f = arxstruc(id, val, nn_f);
% 1st line from v_f -> MSE 
% 2nd line from v_f -> na
% 3rd line from v_f -> nb
% last line from v_f -> delay

%% prediction
na_final = 16; nb_final = 16;
L = length(u_id); % the matrix will have L rows

% identification 
% we initialize the a matrix with zero values
a_pred = zeros(L, 2*na_final); % a_pred = phi matrix id identification

% build the first part of the matrix with y values
for g = 1:L
    l = 1;
    while (l<=na_final && l<=g) 
    a_pred(g,l) = y_id(g - l + 1);
    l = l + 1;
    end
end

% build the second part of the matrix with u values
for g = 1:L
    l = na_final + 1;
    while (l <= 2*na_final+1 && l <= g)
        a_pred(g, l) = u_id(g - l + 1);
        l = l + 1;
    end
end

theta_pred = a_pred\y_id; % the theta used for prediction

% validation
% we initialize the b matrix with zeros values
b_pred = zeros(L, 2*na_final); % b_pred = phi matrix val identification

% build the first part of the matrix with y values
for g = 1:L
    l = 1;
    while (l <= na_final && l <= g) 
    b_pred(g, l) = y_val(g - l + 1);
    l = l + 1;
    end
end

% build the second part of the matrix with u values
for g = 1:L
    l = na_final + 1;
    while (l <= 2*na_final+1 && l <= g)
        b_pred(g, l) = u_val(g - l + 1);
        l = l + 1;
    end
end

% y hat
yh_pred = b_pred * theta_pred;

%% simulation

a_sim = zeros(L, 2*na_final); % a_sim = phi matrix simulation

% build the first part of the matrix with y values
for g = 1:L
    l = 1;
    while (l <= na_final && l <= g) 
        a_sim(g, l) = yh_pred(g - l + 1);
        l = l + 1;
    end
end

% build the second part of the matrix with u values
for g = 1:L
    l = na_final + 1;
    while (l <= 2*na_final+1 && l <= g)
        a_sim(g, l) = u_val(g - l + 1);
        l = l + 1;
    end 
end

yh_sim = a_sim * theta_pred;

%% final plots
figure
plot(yh_pred, 'LineWidth', 2, 'Color', 'yellow'); hold on
plot(yh_sim, 'LineWidth', 1.1); hold on
%plot(y_val, 'Color', 'blue');
legend('y hat prediction', 'y hat simulation'); grid;
title('Final plot with y hat: prediction vs simulation');
xlabel('time'), ylabel('yh');

figure
model_arx = arx(id, [na_final nb_final nk_f]);
compare(model_arx, val);