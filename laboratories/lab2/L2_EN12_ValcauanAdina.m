%% LAB 2 -> Linear Regression for Function Approximation
% Valcauan Adina - Diana 30331/2

close all; 
load('lab2_07.mat')

for n = 2:20

%% id
plot(id.X, id.Y)

% finding phi
N = length(id.X); %n = 5;

for i = 1:N
    for j = 1:n
        phi_id(i,j) = id.X(i)^(j-1);
    end
end

% finding tetha
theta_id = phi_id \ transpose(id.Y);

%% val
plot(val.X, val.Y)

% finding phi
N = length(val.X); %n = 5;

for i = 1:N
    for j = 1:n
        phi_val(i,j) = val.X(i)^(j-1);
    end
end


Y1 = phi_val * theta_id;

%% MSE
s = 0;
for i=1:n
s = s + (val.Y(i)-Y1(i))^2;
end

mse(n) = 1 / n * s;

end

%% plot MSE
n = 2:20;
plot(n,mse(n))
%% 
nMSE = 5;
% id
plot(id.X, id.Y)

% finding phi
N = length(id.X);

for i = 1:N
    for j = 1:nMSE
        phi_id(i,j) = id.X(i)^(j-1);
    end
end

% finding tetha
theta_id = phi_id \ transpose(id.Y);

% val
plot(val.X, val.Y)

% finding phi
N = length(val.X); 

for i = 1:N
    for j = 1:nMSE
        phi_val(i,j) = val.X(i)^(j-1);
    end
end

Y1 = phi_val * theta_id;

% getting the graphic
plot(val.X, val.Y, 'DisplayName','True values'); hold
plot(val.X, Y1, 'DisplayName', 'Approximated');
xlabel('x'); ylabel('y'); legend; title('Final plot')


