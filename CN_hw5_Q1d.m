clear; clc; close all;

eps_values = [0.005; 0.01; 0.02; 0.05; 0.1; 0.2];
tau_values = [2; 3; 5; 7.5; 10; 15; 20];
cases = [10*ones(6,1), eps_values;
    tau_values, 0.02*ones(7,1)];

period = zeros(13,1);
amplitude = zeros(13,1);
lag = zeros(13,1);
f = @(x) 1./(1 + exp(-4*x));
pI = log(0.1/0.9)/4 - 1;
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

for k = 1:13
    tauI = cases(k,1);
    epsilon = cases(k,2);
    wHopf = (1 + 5/tauI)/0.36;
    wEE = (1 + epsilon)*wHopf;

    % Keep the equilibrium at E = I = 0.1.
    pE = log(0.1/0.9)/4 - 0.1*wEE + 1;
    rhs = @(t,y) [(-y(1) + f(wEE*y(1) - 10*y(2) + pE))/5;
        (-y(2) + f(10*y(1) + pI))/tauI];

    % Allow enough time for the smallest epsilon to settle.
    [~,y] = ode45(rhs,[0 20000],[0.101; 0.1],options);
    [t,y] = ode45(rhs,0:0.01:500,y(end,:)',options);
    E = y(:,1);
    I = y(:,2);

    % A peak is higher than its two neighboring samples.
    peaksE = find(E(2:end-1) > E(1:end-2) & ...
        E(2:end-1) >= E(3:end)) + 1;
    peaksI = find(I(2:end-1) > I(1:end-2) & ...
        I(2:end-1) >= I(3:end)) + 1;
    tE = t(peaksE);
    tI = t(peaksI);

    period(k) = mean(diff(tE));
    amplitude(k) = max(E) - min(E);

    % Match each E peak to the next I peak.
    delays = zeros(length(tE)-1,1);
    for j = 1:length(delays)
        delays(j) = tI(find(tI > tE(j),1)) - tE(j);
    end
    lag(k) = 360*mean(delays)/period(k);

    if tauI == 10 && epsilon == 0.1
        trace_t = t(t <= 100);
        trace_y = y(t <= 100,:);
    end
end

fprintf('(i) tau_I = 10 ms\n');
fprintf(' epsilon   Period(ms)    Amplitude     Lag(deg)\n');
fprintf('%8.3g %12.4g %12.4g %12.3g\n', ...
    [eps_values, period(1:6), amplitude(1:6), lag(1:6)]');

fit = polyfit(log(eps_values),log(amplitude(1:6)),1);
fprintf('Fitted slope = %.4g\n',fit(1));

figure;
loglog(eps_values,amplitude(1:6),'o', ...
    eps_values,exp(polyval(fit,log(eps_values))),'--','LineWidth',1.5);
xlabel('\epsilon');
ylabel('Amplitude: max E - min E');
title(sprintf('\\tau_I = 10 ms, slope = %.4f',fit(1)));
legend('Simulation','Fitted line','Location','northwest');
grid on;

f_sim = 1000./period(7:13);
f_linear = 1000/(2*pi)*sqrt(2.592./tau_values - 1./tau_values.^2);

fprintf('\n(ii) epsilon = 0.02\n');
fprintf('tau_I(ms)   Simulated(Hz)   Linear(Hz)\n');
fprintf('%9.3g %15.4g %12.4g\n',[tau_values,f_sim,f_linear]');

figure;
plot(trace_t,trace_y,'LineWidth',1.5);
xlabel('Time after settling (ms)');
ylabel('Activity (dimensionless)');
title('\tau_I = 10 ms, \epsilon = 0.1');
legend('E','I');
grid on;