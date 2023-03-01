%% PROJECT -> Time series modeling using Fourier basis function

% Osiac Maria - Teodora
% Valcauan Adina - Diana
% 30331 / 2

clc; clear
%% loading the data from the .mat file
load('product_30.mat');

figure
plot(time,y, 'blue', 'LineWidth', 1.2); grid
title("The plot for the data")
ylabel('y'), xlabel('t')


% the values we are going to work with
length1 = length(y);
floor_length = floor(length1 * 0.8); % we choose 80% of the data for identification

% the data for identification
yId = y(1:floor_length);
timeId = time(1:floor_length);
lengthId = length(yId);

figure
plot(timeId, yId, 'blue', 'LineWidth', 1.2);
title("The plot for identification data"); grid
ylabel('y'), xlabel('t')


% the data for validation
yVal = y(floor_length:length1);
timeVal = time(floor_length:length1);
lengthVal = length(yVal);

figure
plot(timeVal, yVal, 'blue', 'LineWidth', 1.2);
title("The plot for validation data"); grid
ylabel('y'), xlabel('t')

%% the algorithm for varying m
% we created a for loop, which helps us having a varying m, from 1 to 7
% this also helps us with the MSE
for m = 1:7

%% identification algorithm
P = 12; % the given period (because there are 12 months

for i = 1:length(yId) % we iterate through out identification data
    a(i,1) = 1; a(i,2) = timeId(i); % we populate the first 2 columns with 1 and the desire variable from time array
    for j = 1:m % we start to populate the rest of the columns
        a(i,2*j+1:2*j+2) = [cos(2*pi*j*timeId(i)/P) sin(2*pi*j*timeId(i)/P)];  % we use pairs of columns in order to optain the desire form
    end
end

theta = a \ yId; % the parameters matrix computation

yH_id = a*theta; % now we get y_hat -> estimated values

figure
plot(timeId,yId, 'LineWidth', 1.1); hold on % here we plot the 2 graphs together in
plot (timeId,yH_id, 'LineWidth', 1.1);
title("The plot for every m using identification data")
legend('y identification', 'y_H identification')
ylabel('y, y_H'), xlabel('t')

%% validation algorithm
% very similar with the identification one
P = 12;

for i = 1:length(yVal)
    b(i,1) = 1; b(i,2) = timeVal(i);
    for j = 1:m
        b(i,2*j+1:2*j+2) = [cos(2*pi*j*timeVal(i)/P) sin(2*pi*j*timeVal(i)/P)];
    end
end
yH_val = b*theta;

figure
plot(timeVal,yVal, 'black', 'LineWidth', 1.1); hold
plot (timeVal,yH_val, 'red', 'LineWidth', 1.1);
title("The plot for every m using validation data")
legend('y validation', 'y_H validation')
ylabel('y, y_H'), xlabel('t')

%% MSE
mse_id(m) = (1 / lengthId) .* sum((yId - yH_id).^2); % the MSE for the identification 
mse_val(m) = (1 / lengthVal) .* sum((yVal - yH_val).^2); % the MSE for the validation

end

% plotting the MSE for identification and verification
n = 1:7;
figure
plot(n, mse_id(n), 'green', 'LineWidth', 2); hold on
plot(n, mse_val(n), 'cyan', 'LineWidth', 2); grid
title("The MSE")
legend('identification', 'validation')
ylabel('mse(m)'), xlabel('m')

%% Plotting the best m
% now that we discover the best value of m, we retrace the algorithm in order
% to see the bigger picture

m_best = 3;
P = 12;

% identification part
for i = 1:length(yId)
    a(i,1) = 1; a(i,2) = timeId(i);
    for j = 1:m_best
        a(i,2*j+1:2*j+2) = [cos(2*pi*j*timeId(i)/P) sin(2*pi*j*timeId(i)/P)];
    end
end

theta = a \ yId;

yH_id = a*theta;
figure
plot(timeId,yId, 'magenta', 'LineWidth', 1.2); hold on
plot (timeId,yH_id, 'black', 'LineWidth', 1.1);
title("The plot for the best m using identification data")
legend('y validation', 'y_H validation')
ylabel('y, y_H'), xlabel('t')

% validation part
for i = 1:length(yVal)
    b(i,1) = 1; b(i,2) = timeVal(i);
    for j = 1:m_best
        b(i,2*j+1:2*j+2) = [cos(2*pi*j*timeVal(i)/P) sin(2*pi*j*timeVal(i)/P)];
    end
end
yH_val = b*theta;

figure
plot(timeVal,yVal, 'r', 'LineWidth', 1.1); hold
plot (timeVal,yH_val, 'g', 'LineWidth', 1.1);
title("The plot for the best m using validation data")
legend('y validation', 'y_H validation')
ylabel('y, y_H'), xlabel('t')



