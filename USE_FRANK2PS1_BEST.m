%% Initial Variable and Parameter Setup and Preallocation

% Model Parameters
beta = 0.96;
gamma = 1.3;
r = 0.04;
sigma = 0.04;
rho = 0.9;  % Assumed value

% Income Grid Setup
Y_n = 5;  % Number of income gridpoints
sd = Y_n/2 - 0.5;  % Standard deviations away from mean
Y = linspace(-sd*sigma, sd*sigma, Y_n);  % Assuming normally distributed income around 0

% Asset Grid Setup
a_n = 1000;  % Number of asset gridpoints
a_max = 4 * exp(max(Y));  % Max of assets
A = linspace(-exp(min(Y))/r, a_max, a_n)';  % Asset Grid Discretization

% Transition Probability Matrix (Assumed Simple Case)
P = ones(Y_n) / Y_n;  % Equal transition probabilities

%% Preallocation
c_choice = zeros(a_n, Y_n, a_n);
utility = c_choice;

for ap = 1:a_n
    c(:,:,ap) = (1+r) * repmat(A, 1, Y_n) + exp(repmat(Y, a_n, 1)) - repmat(A(ap), a_n, Y_n);
end
c(c < 0) = 0;  % Adjusting for negative consumption

if gamma == 1
    utility = log(c);
else
    utility(c == 0) = -Inf;
    utility(c > 0) = c(c > 0).^(1-gamma) / (1-gamma);
end

%% VFI Preallocations and Tolerances
tol = 1e-9;
maxits = 1e4;
V0 = zeros(a_n, Y_n);  % Initial Guess of the Value Function
V1 = V0;

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

%% Plots
figure(1)
plot(A, V1(:,1), A, V1(:,round(Y_n/2)), A, V1(:,Y_n))
xlabel('Assets'), ylabel('Value'), title('Value Function')
legend('Minimum Income', 'Medium Income', 'Maximum Income', 'Location', 'SouthOutside', 'Orientation', 'Horizontal')

figure(2)
plot(A, c_policy(:,1), A, c_policy(:,round(Y_n/2)), A, c_policy(:,Y_n))
xlabel('Assets'),
