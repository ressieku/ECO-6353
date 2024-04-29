beta = 0.96;
r = 0.04;
a_n = 50000;
y_n = 5;
rho = 0.92;
sigma = 0.05; % Define standard deviation sigma
sd = 3;
a_max = 10;
sims = 10000;
small_positive = 10^(-4);

% Create transition probability matrix P
P = zeros(y_n);
for i = 1:y_n
    for j = 1:y_n
        P(i, j) = normpdf(j, rho * i, sigma); % Calculate probabilities using normal PDF
    end
    P(i, :) = P(i, :) / sum(P(i, :)); % Normalize rows to sum to 1
end

% Generate income levels y
y = linspace(-sd * sqrt(1 / (1 - rho^2)), sd * sqrt(1 / (1 - rho^2)), y_n);

a_min = -exp(y(1)) / r;

a = linspace(a_min, a_max, a_n)';
[A, Y] = ndgrid(a, y);
w_orig = linspace(a_min + exp(y(1)), a_max + exp(y(end)), a_n)';
w1 = A + Y;

tol = 10^(-9);
maxits = 10^7;

gamma = 3.5; % Define gamma
w = repmat(w_orig, 1, y_n);
c2 = ones(a_n, y_n); % Initialize consumption with ones
c1 = max(w1, c2);
dif = 1;
count = 1;

tic    
while count < maxits && dif > tol
    for ys = 1:y_n
        c2(:, ys) = (beta * (1 + r) * interp1(w_orig, a, w1, 'linear', 'extrap').^(-gamma) * P(ys, :)').^(-1 / gamma); % Use interp1 to interpolate consumption
    end
    c2(c2 < 0) = small_positive; % Ensure non-negative consumption
    w = A / (1 + r) + c2; % Update wealth
    dif = max(max(abs(c1 - c2))); % Calculate difference
    count = count + 1;
    c1 = c2; % Update previous consumption
end
toc

figure(1)
subplot(2, 1, 1)
for i = 1:y_n
    plot(w1(:, i), interp1(real(w(:, i)), real(c1(:, i)), w1(:, i), 'pchip', 'extrap')) % Interpolate consumption for plotting
    hold on
end
hold off
legend('lowest', 'y2', 'steady state', 'y4', 'highest')
xlabel('Cash-on-Hand')
ylabel('Consumption')

subplot(2, 1, 2)
plot(w1(:, 1), interp1(real(w(:, 1)), real(c1(:, 1)), w1(:, 1), 'pchip', 'extrap')) % Interpolate consumption for plotting
legend('lowest income')
xlabel('Cash-on-Hand')
ylabel('Consumption')

figure(2)
subplot(2, 1, 1)
for i = 1:y_n
    plot(a, interp1(real(w(:, i)), real(c1(:, i)), a, 'pchip', 'extrap')) % Interpolate consumption for plotting
    hold on
end
hold off
legend('lowest', 'y2', 'steady state', 'y4', 'highest')
xlabel('Assets')
ylabel('Consumption')

subplot(2, 1, 2)
plot(a, interp1(real(w(:, 1)), real(c1(:, 1)), a, 'pchip', 'extrap')) % Interpolate consumption for plotting
legend('lowest income')
xlabel('Assets')
ylabel('Consumption')

y_sim = simulate(dtmc(P), sims);
w_sim = zeros(1, sims);

for t = 1:sims
    c_sim(t) = interp1(real(w(:, y_sim(t))), real(c1(:, y_sim(t))), real(w_sim(t)), 'spline', min(min(c1))); % Interpolate consumption for simulation
    if t < sims
        w_sim(t + 1) = (1 + r) * (w_sim(t) - c_sim(t)) + exp(y(y_sim(t + 1))); % Update wealth for simulation
    end
end

figure(3)
subplot(3, 1, 1)
plot(sims / 2 + 1:sims, exp(y(y_sim(sims / 2 + 1:sims))), sims / 2 + 1:sims, ones(sims / 2, 1))
xlabel('Time')
ylabel('Income')

subplot(3, 1, 2)
plot(sims / 2 + 1:sims, c_sim(sims / 2 + 1:sims))
xlabel('Time')
ylabel('Consumption')

subplot(3, 1, 3)
plot(sims / 2 + 1:sims, w_sim(sims / 2 + 1:sims))
xlabel('Time')
ylabel('Cash-On-Hand')


% Part (b) No Borrowing

beta = 0.96;
r = 0.04;
a_n = 50000;
y_n = 5;
rho = 0.92;
sigma = 0.05; % Define standard deviation sigma
sd = 3;
a_max = 10;
sims = 10000;
small_positive = 10^(-4);

P = zeros(y_n);
for i = 1:y_n
    for j = 1:y_n
        P(i, j) = normpdf(j, rho * i, sigma); % Calculate probabilities using normal PDF
    end
    P(i, :) = P(i, :) / sum(P(i, :)); % Normalize rows to sum to 1
end

y = linspace(-sd * sqrt(1 / (1 - rho^2)), sd * sqrt(1 / (1 - rho^2)), y_n);

a_min = 0; % Minimum asset level

a = linspace(a_min, a_max, a_n)';
[A, Y] = ndgrid(a, y);
w_orig = linspace(a_min + exp(y(1)), a_max + exp(y(end)), a_n)';
w1 = A + Y;

tol = 10^(-9);
maxits = 10^7;

gamma = 3.5; % Define gamma
w = repmat(w_orig, 1, y_n);
c2 = ones(a_n, y_n); % Initialize consumption with ones
c1 = max(w1, c2);
dif = 1;
count = 1;

tic    
while count < maxits && dif > tol
    for ys = 1:y_n
        c2(:, ys) = (beta * (1 + r) * interp1(w_orig, a, w1, 'linear', 'extrap').^(-gamma) * P(ys, :)').^(-1 / gamma); % Use interp1 to interpolate consumption
    end
    c2(c2 < 0) = small_positive; % Ensure non-negative consumption
    w = A / (1 + r) + c2; % Update wealth
    dif = max(max(abs(c1 - c2))); % Calculate difference
    count = count + 1;
    c1 = c2; % Update previous consumption
end
toc

figure(4)
subplot(5, 1, 1)
for i = 1:y_n
    plot(w1(:, i), interp1(real(w(:, i)), real(c1(:, i)), w1(:, i), 'pchip', 'extrap')) % Interpolate consumption for plotting
    hold on
end
hold off
legend('lowest', 'y2', 'steady state', 'y4', 'highest')
xlabel('Cash-on-Hand')
ylabel('Consumption')

subplot(5, 1, 2)
plot(w1(:, 1), interp1(real(w(:, 1)), real(c1(:, 1)), w1(:, 1), 'pchip', 'extrap')) % Interpolate consumption for plotting
legend('lowest income')
xlabel('Cash-on-Hand')
ylabel('Consumption')

figure(5)
subplot(5, 1, 1)
for i = 1:y_n
    plot(a, interp1(real(w(:, i)), real(c1(:, i)), a, 'pchip', 'extrap')) % Interpolate consumption for plotting
    hold on
end
hold off
legend('lowest', 'y2', 'steady state', 'y4', 'highest')
xlabel('Assets')
ylabel('Consumption')

subplot(5, 1, 2)
plot(a, interp1(real(w(:, 1)), real(c1(:, 1)), a, 'pchip', 'extrap')) % Interpolate consumption for plotting
legend('lowest income')
xlabel('Assets')
ylabel('Consumption')

y_sim = simulate(dtmc(P), sims);
w_sim = zeros(1, sims);

for t = 1:sims
    c_sim(t) = interp1(real(w(:, y_sim(t))), real(c1(:, y_sim(t))), real(w_sim(t)), 'spline', min(min(c1))); % Interpolate consumption for simulation
    if t < sims
        w_sim(t + 1) = (1 + r) * (w_sim(t) - c_sim(t)) + exp(y(y_sim(t + 1))); % Update wealth for simulation
    end
end

figure(6)
subplot(6, 1, 1)
plot(sims / 2 + 1:sims, exp(y(y_sim(sims / 2 + 1:sims))), sims / 2 + 1:sims, ones(sims / 2, 1))
xlabel('Time')
ylabel('Income')

subplot(6, 1, 2)
plot(sims / 2 + 1:sims, c_sim(sims / 2 + 1:sims))
xlabel('Time')
ylabel('Consumption')

subplot(6, 1, 3)
plot(sims / 2 + 1:sims, w_sim(sims / 2 + 1:sims))
xlabel('Time')
ylabel('Cash-On-Hand')


% Part (c) Sigma doubled -> 0.1

beta = 0.96;
r = 0.04;
a_n = 50000;
y_n = 5;
rho = 0.92;
sigma = 0.1; % Define standard deviation sigma
sd = 3;
a_max = 10;
sims = 10000;
small_positive = 10^(-4);

P = zeros(y_n);
for i = 1:y_n
    for j = 1:y_n
        P(i, j) = normpdf(j, rho * i, sigma); % Calculate probabilities using normal PDF
    end
    P(i, :) = P(i, :) / sum(P(i, :)); % Normalize rows to sum to 1
end

y = linspace(-sd * sqrt(1 / (1 - rho^2)), sd * sqrt(1 / (1 - rho^2)), y_n);

a_min = -exp(y(1)) / r;

a = linspace(a_min, a_max, a_n)';
[A, Y] = ndgrid(a, y);
w_orig = linspace(a_min + exp(y(1)), a_max + exp(y(end)), a_n)';
w1 = A + Y;

tol = 10^(-9);
maxits = 10^7;

gamma = 3.5; % Define gamma
w = repmat(w_orig, 1, y_n);
c2 = ones(a_n, y_n); % Initialize consumption with ones
c1 = max(w1, c2);
dif = 1;
count = 1;

tic    
while count < maxits && dif > tol
    for ys = 1:y_n
        c2(:, ys) = (beta * (1 + r) * interp1(w_orig, a, w1, 'linear', 'extrap').^(-gamma) * P(ys, :)').^(-1 / gamma); % Use interp1 to interpolate consumption
    end
    c2(c2 < 0) = small_positive; % Ensure non-negative consumption
    w = A / (1 + r) + c2; % Update wealth
    dif = max(max(abs(c1 - c2))); % Calculate difference
    count = count + 1;
    c1 = c2; % Update previous consumption
end
toc

figure(7)
subplot(8, 1, 1)
for i = 1:y_n
    plot(w1(:, i), interp1(real(w(:, i)), real(c1(:, i)), w1(:, i), 'pchip', 'extrap')) % Interpolate consumption for plotting
    hold on
end
hold off
legend('lowest', 'y2', 'steady state', 'y4', 'highest')
xlabel('Cash-on-Hand')
ylabel('Consumption')

subplot(8, 1, 2)
plot(w1(:, 1), interp1(real(w(:, 1)), real(c1(:, 1)), w1(:, 1), 'pchip', 'extrap')) % Interpolate consumption for plotting
legend('lowest income')
xlabel('Cash-on-Hand')
ylabel('Consumption')

figure(8)
subplot(8, 1, 1)
for i = 1:y_n
    plot(a, interp1(real(w(:, i)), real(c1(:, i)), a, 'pchip', 'extrap')) % Interpolate consumption for plotting
    hold on
end
hold off
legend('lowest', 'y2', 'steady state', 'y4', 'highest')
xlabel('Assets')
ylabel('Consumption')

subplot(8, 1, 2)
plot(a, interp1(real(w(:, 1)), real(c1(:, 1)), a, 'pchip', 'extrap')) % Interpolate consumption for plotting
legend('lowest income')
xlabel('Assets')
ylabel('Consumption')

y_sim = simulate(dtmc(P), sims);
w_sim = zeros(1, sims);

for t = 1:sims
    c_sim(t) = interp1(real(w(:, y_sim(t))), real(c1(:, y_sim(t))), real(w_sim(t)), 'spline', min(min(c1))); % Interpolate consumption for simulation
    if t < sims
        w_sim(t + 1) = (1 + r) * (w_sim(t) - c_sim(t)) + exp(y(y_sim(t + 1))); % Update wealth for simulation
    end
end

figure(9)
subplot(9, 1, 1)
plot(sims / 2 + 1:sims, exp(y(y_sim(sims / 2 + 1:sims))), sims / 2 + 1:sims, ones(sims / 2, 1))
xlabel('Time')
ylabel('Income')

subplot(9, 1, 2)
plot(sims / 2 + 1:sims, c_sim(sims / 2 + 1:sims))
xlabel('Time')
ylabel('Consumption')

subplot(9, 1, 3)
plot(sims / 2 + 1:sims, w_sim(sims / 2 + 1:sims))
xlabel('Time')
ylabel('Cash-On-Hand')


% Part (d) Correlogram

[clgm, lags] = xcorr(y_sim, c_sim, 4);

% Part (e) Plotting the Correlogram

figure(10)
plot(lags, clgm);
xlabel("lags")
ylabel("Correlation Vector")
title("Correlogram between Simulated Income and Consumption Series")
legend("Correlation at 4 lags")
