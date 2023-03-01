%% LAB 7 -> Pseudo-random binary sequences
% Valcauan Adina - Diana 30331/2
close all;

%% the algorithm for m = 3
clc; clear; 

% loaded the data
load('uval.mat');
length_u = length(u);

m = 3; a = 0.5; b = 1; % given in the lab pdf
u_prbs = PRBS(length_u, m, a, b); % calling the created function
% u_prbs = idinput(length_u, 'prbs', [], [a, b]); % used for verification purposes

% the index is the same like in the last lab: 1
id_data = system_simulator(1, u_prbs); % identification
val_data = system_simulator(1, u); % validation

% plotting the resulted data
figure
plot(id_data, 'magenta')
title("Identification data for m = 3", "Color", 'blue')
figure
plot(val_data, 'magenta')
title("Validation data for m = 3", "Color", 'blue')

% doing the simulated resonse comparison
na = 15; nb = 15; nk = 1;
figure
model_arx = arx(id_data, [na nb nk]);
compare(model_arx, val_data);
title('Simulated Response Comparison for m = 3')

% the period for the maximum-length PRBS
period_m3 = 2^(m)-1;

%% the algorithm for m = 10
clc; clear;

% loaded the data
load('uval.mat');
length_u = length(u);

m = 10; a = 0.5; b = 1; % given in the lab pdf
u_prbs = PRBS(length_u, m, a, b); % calling the function for a different m
% u_prbs = idinput(length_u, 'prbs', [], [a, b]); % used for verification purposes

% the index is the same like in the last lab: 1
id_data = system_simulator(1, u_prbs); % identification
val_data = system_simulator(1, u); % validation

% plotting the resulted data
figure
plot(id_data, 'red')
title("Identification data for m = 10", "Color", 'blue')
figure
plot(val_data, 'red')
title("Validation data for m = 10", "Color", 'blue')

% doing the simulated resonse comparison
na = 15; nb = 15; nk = 1;
figure
model_arx = arx(id_data, [na nb nk]);
compare(model_arx, val_data);
title('Simulated Response Comparison for m = 10')

% the period for the maximum-length PRBS
period_m10 = 2^(m)-1;
%% the PRBS function needed at step no 2
function [u_prbs] = PRBS(length_u, m, a, b)

% building the A matrix from the course
a_array = zeros(1, m); % the array represents the first line from the matrix
% the switch is used for the correct coefficients
switch m
    case 3
        a_array(1) = 1; a_array(3) = 1;
    case 4
        a_array(1) = 1; a_array(4) = 1;
    case 5
        a_array(2) = 1; a_array(5) = 1;
    case 6 
        a_array(1) = 1; a_array(6) = 1;
    case 7
        a_array(1) = 1; a_array(7) = 1;
    case 8
        a_array(1) = 1; a_array(2) = 1; a_array(7) = 1; a_array(8) = 1;
    case 9
        a_array(4) = 1; a_array(9) = 1;
    case 10
        a_array(3) = 1; a_array(10) = 1;  
end

% the A matrix
matrix(1, :) = a_array(1:m);
matrix(2:m, 1:m-1) = eye(m-1);
matrix(2:m, m) = zeros(m-1, 1);

% the C matrix from the course
c_matrix = zeros(1, m-1);
c_matrix(m) = 1;

u_prbs_unshifted = zeros(length_u, 1);
u_prbs = zeros(length_u, 1);
% trying to do the state space representation
for i = 1:length_u
    X_k = ones(m, 1); % x(k+0)
    for j = 1:i-1
    X_k_plus_1 = mod(matrix*X_k, 2); % x(k+1)
    X_k = X_k_plus_1;
    end
    u_prbs_unshifted(i) = c_matrix * X_k; % the u_prbs with 0 and 1
    u_prbs(i) = a + (b - a) * u_prbs_unshifted(i); % shifting
end

end