clear; clc; close all;
rng(1);

mu = 0.2;
sigma = 0.5;
I0 = 5;
dt = 0.002;
N = 10000;

theta = 0.5:0.5:4;
lambda = [0, 0.2];

accuracy = zeros(2, length(theta));
meanTime = zeros(size(accuracy));

%% Simulate the two-population model
for j = 1:length(lambda)
    ws = (1 + lambda(j))/2;
    wc = ws;

    for k = 1:length(theta)
        r1 = I0*ones(N,1);
        r2 = r1;

        t = 0;
        correct = 0;
        totalTime = 0;

        while ~isempty(r1)
            n = length(r1);

            % Both updates use the OLD rates and independent noise.
            dr1 = ((ws-1)*r1 - wc*r2 + I0 + mu/2)*dt + sigma*sqrt(dt)*randn(n,1);

            dr2 = ((ws-1)*r2 - wc*r1 + I0 - mu/2)*dt + sigma*sqrt(dt)*randn(n,1);
                  

            r1 = r1 + dr1;
            r2 = r2 + dr2;
            t = t + dt;

            d = r1 - r2;
            finished = abs(d) >= theta(k);

            correct = correct + sum(d >= theta(k));
            totalTime = totalTime + t*sum(finished);

            % Remove trials as soon as either threshold is reached.
            r1(finished) = [];
            r2(finished) = [];
        end

        accuracy(j,k) = correct/N;
        meanTime(j,k) = totalTime/N;
    end
end

%% Predictions from parts (b) and (c)
h = linspace(theta(1), theta(end), 200);

% Perfect integrator: accuracy and mean decision time
P0 = 1./(1 + exp(-mu*h/sigma^2));
T0 = (h/mu).*tanh(mu*h/(2*sigma^2));

% Unstable circuit: part (c) gives an accuracy prediction
g = @(y) exp(-(lambda(2)*y.^2 + 2*mu*y)/(2*sigma^2));
P02 = zeros(size(h));

for k = 1:length(h)
    P02(k) = integral(g, -h(k), 0) ...
             / integral(g, -h(k), h(k));
end

%% Plot accuracy
figure;
plot(h, P0, '-', h, P02, '--', 'LineWidth', 1.5);
hold on;
plot(theta, accuracy(1,:), 'o', theta, accuracy(2,:), 's');
grid on;

xlabel('Threshold \theta');
ylabel('Probability of a correct choice');
title('\mu=0.2, \sigma=0.5, I_0=5, dt=0.002, N=10^4');
legend('\lambda=0: theory', '\lambda=0.2: theory', ...
       '\lambda=0: simulation', '\lambda=0.2: simulation', ...
       'Location', 'southeast');

%% Plot mean decision time
figure;
plot(h, T0, '-', theta, meanTime(1,:), 'o', ...
     theta, meanTime(2,:), 's-', 'LineWidth', 1.5);
grid on;

xlabel('Threshold \theta');
ylabel('Mean decision time (rate time constants)');
title('\mu=0.2, \sigma=0.5, I_0=5, dt=0.002, N=10^4');
legend('\lambda=0: theory', '\lambda=0: simulation', ...
       '\lambda=0.2: simulation', 'Location', 'northwest');

%% Compare at equal predicted accuracy
% Accuracy for lambda=0.2, theta=4
pMatch = P02(end);

% Invert the perfect-integrator accuracy formula
thetaMatch = (sigma^2/mu)*log(pMatch/(1-pMatch));
Tmatch = (thetaMatch/mu)*tanh(mu*thetaMatch/(2*sigma^2));

fprintf('\nUnstable circuit: lambda=0.2, theta=4\n');
fprintf('Accuracy: %.4f (theory), %.4f (simulation)\n', ...
        pMatch, accuracy(2,end));
fprintf('Mean decision time: %.4g (simulation)\n', meanTime(2,end));

fprintf('\nPerfect integrator at the same predicted accuracy:\n');
fprintf('Threshold: %.4g\n', thetaMatch);
fprintf('Accuracy: %.4f\n', pMatch);
fprintf('Mean decision time: %.4g (theory)\n', Tmatch);