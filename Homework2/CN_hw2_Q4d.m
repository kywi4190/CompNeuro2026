clear; clc; close all;

% time constant (s)
tau = 0.02;
% mean voltage (mV)
mu = 10;
nuTau = [2 20 200];

% time step (s)
dt = 0.0001;
% simulation time (s)
T = 100;
N = round(T/dt);

results = zeros(3,7);

for j = 1:3

    % Choose nu and w so that mean voltage stays at 10 mV
    nu = nuTau(j)/tau;
    w = mu/nuTau(j);

    % Number of Poisson events in each time step
    events = poissrnd(nu*dt, N, 1);

    % Simulate voltage
    V = zeros(N,1);
    decay = exp(-dt/tau);

    for k = 2:N
        V(k) = decay*V(k-1) + w*events(k);
    end

    % Ignore the first second while V settles
    V = V(round(1/dt):end);

    % Measured statistics
    m = mean(V);
    v = var(V,1);
    s = mean((V-m).^3)/v^(3/2);

    % Theoretical statistics
    vTheory = nu*w^2*tau/2;
    sTheory = 2*sqrt(2)/(3*sqrt(nuTau(j)));

    results(j,:) = [nuTau(j), m, mu, v, vTheory, s, sTheory];

    % Histograms for nu*tau = 2 and 200
    if nuTau(j) == 2 || nuTau(j) == 200
        figure
        histogram(V,60,'Normalization','pdf')
        hold on

        x = linspace(min(V),max(V),300);
        G = 1/sqrt(2*pi*vTheory) * ...
            exp(-(x-mu).^2/(2*vTheory));

        plot(x,G,'LineWidth',2)

        xlabel('V (mV)')
        ylabel('Probability density')
        title(sprintf('\\nu\\tau = %g',nuTau(j)))
        legend('Shot noise','Gaussian')
    end
end

% Display the measured values against theoretic values
table(results(:,1), results(:,2), results(:,3), ...
    results(:,4), results(:,5), results(:,6), results(:,7), ...
    'VariableNames', {'nuTau','Mean','MeanTheory', ...
    'Variance','VarianceTheory','Skewness','SkewTheory'})
