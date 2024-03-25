%% Initial Variable and Parameter Setup and Preallocation

% Model Parameters
beta = 0.96;
gamma = 1.3;
r = 0.04;
sigma = 0.04;
rho = 0.9;  % Assumed value

% Income Grid Setup
Y_n = 5;
sd = Y_n / 2 - 0.5; 
Y = linspace(-sd * sigma, sd * sigma, Y_n); 

% Asset Grid Setup
a_n = 1000;
a_max = 4 * exp(max(Y)); 
A = linspace(-exp(min(Y)) / r, a_max, a_n)';  

% Transition Probability Matrix
P = ones(Y_n) / Y_n;  

% Preallocation
c_choice = zeros(a_n, Y_n, a_n);
utility = c_choice;

%% Main VFI Loop
count = 0;
dif = inf;
while dif > tol && count < maxits
    for y = 1:Y_n
        for ap = 1:a_n
            V_candidate = utility(:,y,ap) + beta * V0 * P(y,:)';
            [V1(:,y), ~] = max(V_candidate, [], 2);
        end
    end
    dif = max(abs(V0(:) - V1(:)));
    V0 = V1;
    count = count + 1;
end

%% Recovery of Consumption Policy Function
% Assuming the a_prime index is already calculated and available
a_prime = ones(size(A));  % Placeholder for asset choice index
c_policy = (1+r) * A + exp(Y(a_prime)) - A(a_prime);

%% Simulating the income process
sims = 1000;  % Number of simulation periods
y_sim = zeros(sims, 1);
y_sim(1) = Y(round(Y_n / 2)); % Start from the middle of the grid

% Simulate income using AR(1) process
for t = 2:sims
    epsilon = sigma * randn; % Random shock with standard deviation sigma
    y_sim(t) = rho * y_sim(t - 1) + epsilon;
end

% Simulate consumption based on simulated income and policy function
% Placeholder for asset choice index, replace with actual data from VFI
a_prime_sim = ones(sims, 1); 
c_sim = (1 + r) * A(a_prime_sim) + exp(y_sim) - A(a_prime_sim);

%% Plots
% Existing plot code...

%% Correlogram calculation
[acor, lags] = xcorr(y_sim, c_sim, 4, 'coeff');  % Calculate correlogram

% Correlogram plot
figure;
% The corrplot function is not standard in MATLAB, you might need to use another plotting function or custom script
% For illustration, here's how you might plot with MATLAB's built-in function
bar(lags, acor);
xlabel('Lag'), ylabel('Cross-correlation'), title('Correlogram between Income and Consumption');
