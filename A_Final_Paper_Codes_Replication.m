%%%%% Replication by Richmond Essieku %%%%%
% Cooper, Russell, John Haltiwanger, and Laura Power. 
% "Machine replacement and the business cycle: lumps and bumps." 
% American Economic Review 89, no. 4 (1999): 921-946.

%%%%%% The Value Function %%%%%%
% The decision model for machine replacement which determines the best action (replace or not)
% based on maximizing the value function V(k, A, e) that has two components:
% V^n(k, A, e) = A*e*f(k) + beta*EV(p*k, A', e') (continuing with existing capital)
% V^r(k, A, e) = A*e*f(k)*lambda - F + beta*EV(1, A', e') (replacing the capital)

% Setting parameter values based on the given problem
AH = 1.25;  % High productivity
AL = 0.75;  % Low productivity
A = [AH, AL];  % Array of productivity levels
delta = 0.1;  % Depreciation rate
beta = 0.9;  % Discount factor
mu = 1;  % Unused parameter in this snippet
rho = 1 - delta;  % Capital retention rate
F = 0.2;  % Cost of replacement
lambda = 0.75;  % Efficiency gain from new capital
N = 20;  % Number of grid points for epsilon
R = 6;  % Replacement period in years
% Transition probabilities for productivity states
P = [0.9, 0.1; 0.1, 0.9];

% Defining the grid for productivity shocks
e_min = 0.4;
e_max = 1.6;
e = linspace(e_min, e_max, N)';  % Epsilon values evenly spaced

% Defining the capital stock array
k = [1, cumprod(rho * ones(1, R - 1))];  % Array of depreciated capital values


% Part 1) Let initialize the following object

% (a) Constructing the three dimensional array following the directive of
% the question.
s = e * k;  % Compute potential revenues for each combination of shock and capital
rev = s * A(1);  % Revenue at high productivity
rev(:,:,2) = s * A(2);  % Revenue at low productivity

% For part(b), part(c), and part (d) will help to initialize the three
% dimensional array  for value functions, which is as follows
V_R = ones(size(rev));  % Initial guess for the value of replacing
V_N = zeros(size(rev));  % Initial guess for the value of not replacing
V0 = max(V_R, V_N);  % The overall initial guess for the value function

tol = 1e-9;  % Convergence tolerance for the iteration

% Part 2) The Value Function Iteration Loop

dif = 1;  % Initialize the difference for the while loop
count = 0;  % Iteration counter
maxits = 1e9;  % Maximum number of iterations

tic  % Start timer
while dif > tol && count < maxits
    z0 = zeros(size(rev));  % Initialize the policy function (replacement decision)
    for i = 1:2  % Loop over each productivity state
        V_R(:,:,i) = rev(:,:,i) * lambda - F + beta * (P(i,1) * repmat(V0(:,1,1), 1, R) + P(i,2) * repmat(V0(:,1,2), 1, R));
        for j = 1:R-1
            V_N(:,j,i) = rev(:,j,i) + beta * (P(i,1) * V0(:,j+1,1) + P(i,2) * V0(:,j+1,2));
        end
    end
    V1 = max(V_R, V_N);  % Update the overall value function
    z0(V1 == V_R) = 1;  % Update the policy (replacement decision)
    dif = max(max(max(abs(V1 - V0))));  % Calculate the maximum difference for convergence
    count = count + 1;  % Increment the counter
    V0 = V1;  % Update value function for the next iteration
end
toc  % Stop timer

% Part 3) Examining the policy function by Updating the cutoff level

R_new = 7;  % New cutoff level after extending the period

k1 = [1, cumprod(rho * ones(1, R_new - 1))];  % New capital stock array

% (a) Recreating the three dimensional array revenue function matrix for the new grid points
s1 = e * k1;  % Potential revenues with new capital stock
rev1 = s1 * A(1);  % High productivity revenue
rev1(:,:,2) = s1 * A(2);  % Low productivity revenue

% (b), (c), and (d) Reinitializing 3D arrays for value functions with new capital stock
V_R1 = ones(size(rev1));  % Initial guess for the value of replacing
V_N1 = zeros(size(rev1));  % Initial guess for the value of not replacing
V00 = max(V_R1, V_N1);  % Overall initial guess for the value function

dif = 1;  % Reset difference for loop
count = 0;  % Reset counter
maxits = 1e9;  % Reset maximum iterations

tic  % Start timer
while dif > tol && count < maxits
    z01 = zeros(size(rev1));  % Initialize policy function for the new capital stock
    for i = 1:2  % Loop over each productivity state for new capital stock
        V_R1(:,:,i) = rev1(:,:,i) * lambda - F + beta * (P(i,1) * repmat(V00(:,1,1), 1, R_new) + P(i,2) * repmat(V00(:,1,2), 1, R_new));
        for j = 1:R_new-1
            V_N1(:,j,i) = rev1(:,j,i) + beta * (P(i,1) * V00(:,j+1,1) + P(i,2) * V00(:,j+1,2));
        end
    end
    V11 = max(V_R1, V_N1);  % Update the overall value function with new capital stock
    z01(V11 == V_R1) = 1;  % Update policy function
    dif = max(max(max(abs(V11 - V00))));  % Calculate maximum difference for convergence
    count = count + 1;  % Increment counter
    V00 = V11;  % Update value function for the next iteration
end
toc  % Stop timer


% Part 4) Ploting the policy function using Spy Plot
% Visualization of the replacement decision over time and productivity states

figure(1)
subplot(2,1,1)
spy(z0(:,:,1)', 'black')  % Spy plot for low aggregate productivity
xlabel('Firm Level productivity')
ylabel('Time series since last replaced')
title('Low Aggregate Productivity Replacement')

subplot(2,1,2)
spy(z0(:,:,2)', 'red')  % Spy plot for high aggregate productivity
xlabel('Firm level  productivity')
ylabel('Time series since last replaced')
title('High Aggregate Productivity Replacement')

% Interpretation of policy function: As capital ages, reaching the threshold (R=6 years),
% firms are more likely to replace old capital especially in higher productivity conditions,
% reflecting a strategic response to balance between capital costs and potential revenue gains.


% Part 5) Ploting the Capital Replacement Hazard Function H(k,A)
% Calculate the probability of replacing capital as it ages

H = zeros(R, 2);  % Initialize hazard function array

for i = 1:2  % Loop through each productivity state
    for j = 1:R  % Loop through each year up to replacement threshold
        H(j, i) = sum(z0(:, j, i)) / N;  % Calculate probability of replacement
    end
end

time = 1:R;  % Time array for plotting

figure(2)
plot(time, H(:, 1), 'red', time, H(:, 2), 'black')  % Plot hazard functions for both productivity states
title('Hazard Function')
xlabel('Time Since Last Replacement')
ylabel('Probability of Replacement')
legend('Low state', 'High State','Location', 'northwest')
xlim([1 R])
ylim([0 1.05])

% Part 6) Time Series Analysis through simulation
% Simulate the impact of capital replacement decisions on firm's output over time

ts = 40;  % Number of time steps for simulation
T = 160;  % Total time for aggregate state simulation
Time = 1:ts;  % Time array for plotting
E = randi(N, 1, ts);  % Randomly generate productivity shocks for each time step
esim = e(E);  % Corresponding productivity values

AT = zeros(1, T);  % Initialize array for aggregate states
AT(1) = 2;  % Start from high productivity state

for i = 2:T  % Markov process to simulate transitions between productivity states
    if AT(i-1) == 1
        AT(i) = randsample(1:2, 1, true, P(1, :));  % Transition from low to either state
    else
        AT(i) = randsample(1:2, 1, true, P(2, :));  % Transition from high to either state
    end
end

Asim = A(AT);  % Corresponding aggregate productivity states

Y = zeros(1, ts);  % Initialize output array
K = zeros(1, ts + 1);  % Initialize capital array
K(1) = 1;  % Start with new capital
sim_K = zeros(1, ts + 1);  % Initialize simulated capital array
sim_K(1) = 1;  % Start with new capital

for i = 1:ts  % Loop over each time step
    simoutput(i) = K(K(i)) * Asim(i) * esim(i) * (1 - z0(E(i), K(i), AT(i))) + z0(E(i), K(i), AT(i)) * (lambda * K(K(i)) * Asim(i) * esim(i) - F);
    if z0(E(i), K(i), AT(i)) == 1
        K(i + 1) = 1;  % Reset capital if replaced
        sim_K(i + 1) = K(K(i + 1));  % Update simulated capital
    else
        K(i + 1) = K(i) + 1;  % Age the capital if not replaced
        sim_K(i + 1) = K(K(i + 1));  % Update simulated capital
    end
end

figure(4)
plot(Time, sim_K(1:ts), 'red--', Time, simoutput, 'black-')  % Plot simulated capital and output
title('Aggregate Investment Fluctuations "Simulated Capital"')
xlabel('Period')
ylabel('Investment Rate')
ylim([0 max(simoutput) + 0.1])
legend('Capital', 'Output', 'location', 'northwest')

% Part 7) Rate of Investment
% Analyze how the investment rate converges over time under fixed aggregate productivity

tt = 50;  % Total time steps for investment rate analysis
n = 6;  % Number of firms
ir = zeros(tt, 1);  % Initialize investment rate array
w = zeros(N, tt + 1);  % Initialize array for number of firms at each capital vintage
w(:, 1) = ones(N, 1);  % Start with one firm at each vintage
A_index = 1;  % Fix productivity state for simplicity

for i = 1:tt  % Loop over each time step
    for j = 1:n  % Loop over each firm
        ir(i) = ir(i) + H(j, A_index) * w(j, i);  % Calculate investment rate
    end
    w(1, i + 1) = ir(i);  % Update number of firms at new capital
    for j = 2:n
        w(j, i + 1) = (1 - H(j - 1, A_index)) * w(j - 1, i);  % Update number of firms at aged capital
    end
end

ir = ir / n;  % Normalize investment rate by number of firms

figure(3)
plot(1:tt, ir, 'red')  % Plot investment rate over time
title('Convergence without Aggregate Shocks "Baseline Parameters"')
xlabel('Period')
ylabel('Investment Rate')

% Part 8) The Markov process
% Simulate investment rate under a Markov process for aggregate state transitions

ir2 = zeros(T - 1, 1);  % Initialize investment rate array for Markov process
w1 = zeros(n, T);  % Initialize array for number of firms at each capital vintage under Markov process
w1(:, 1) = ones(n, 1);  % Start with one firm at each vintage

for i = 1:T - 1  % Loop over each time step
    for j = 1:n  % Loop over each firm
        ir2(i) = ir2(i) + H(j, AT(i)) * w1(j, i);  % Calculate investment rate under Markov process
    end
    w1(1, i + 1) = ir2(i);  % Update number of firms at new capital
    for j = 2:n
        w1(j, i + 1) = (1 - H(j - 1, AT(i))) * w1(j - 1, i);  % Update number of firms at aged capital
    end
end

ir2 = ir2 / n;  % Normalize investment rate by number of firms

figure('Name', '3 and 4');
[y, line1, line2] = plotyy(1:T - 1, ir2, 1:T, Asim);  % Plot investment rate and aggregate state
line1.Color = 'red';  % Set color for investment rate line
line2.Color = 'black';  % Set color for aggregate state line
ylabel(y(1), 'Investment Rate', 'Color', 'black')  % Label for investment rate
ylabel(y(2), 'Aggregate State', 'Color', 'black')  % Label for aggregate state
xlabel('Period', 'FontWeight', 'bold')  % Bold label for period
line2.LineStyle = '-.';  % Dotted line for aggregate state
line2.Marker = '*';  % Marker for aggregate state
title('Aggregate Investment Fluctuations - Baseline Simulation')

% Part 9) Please see the implication described in the attached PDF file

