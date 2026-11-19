clear; clc; close all;

rng(1);

pH = 0.5;
rho = 0.5;
tau = 5;
N = 4e5;

R = @(T) pH*(1-exp(-rho*T)) ./ (pH*(1-exp(-rho*T))/rho + (1-pH)*T + tau);
lambda = @(T) rho*pH*exp(-rho*T) ./ (pH*exp(-rho*T) + 1-pH);
Rknown = pH/(pH/rho + tau);

Tstar = fzero(@(T) R(T)-lambda(T), [0, 10]);
Rstar = R(Tstar);

% (i) Simulate 400,000 visits for each waiting time.
fprintf('   T        Measured      R(T)\n');
for T = [2, Tstar, 10]
    high = rand(N,1) < pH;
    % Exponential waiting times
    findTime = -log(rand(N,1))/rho;
    % Low patches have no food
    findTime(~high) = Inf;

    food = findTime <= T;
    visitTime = min(findTime,T) + tau;
    measured = sum(food)/sum(visitTime);

    fprintf('%7.4f     %.6f     %.6f\n', T, measured, R(T));
end

% (ii) Plot the rates and mark the optimum.
Tplot = linspace(0,15,500);
figure;
plot(Tplot, R(Tplot), 'b', 'LineWidth', 1.5);
hold on;
yline(Rknown, 'k--', 'LineWidth', 1.5);
plot(Tplot, lambda(Tplot), 'r', 'LineWidth', 1.5);
plot(Tstar, Rstar, 'ko', 'MarkerFaceColor', 'k');
text(Tstar+0.3, Rstar, sprintf('T^* = %.3f, R^* = %.5f', Tstar, Rstar), 'VerticalAlignment', 'bottom');
xlabel('Waiting time T');
ylabel('Food per unit time');
legend('$R(T)$', '$R_{\mathrm{known}}$', '$\tilde{\lambda}(T)$', '$(T^*,R^*)$', 'Interpreter', 'latex');
grid on;