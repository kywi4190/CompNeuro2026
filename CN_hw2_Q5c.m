clear; clc; close all;
rng(1);

mu = 1;
theta = 1;
sigma2_list = [0.25 1 4];

dt = 2e-4;
T = 3000;
Tdiscard = 50;

N = round(T/dt);
Ndiscard = round(Tdiscard/dt);

figure;

for k = 1:length(sigma2_list)

    sigma2 = sigma2_list(k);

    u = 0;
    spikes = 0;
    below = 0;
    count = 0;

    % Only save every 20th point for the histogram
    samples = [];
    sampleEvery = 20;

    for n = 1:N

        % Euler-Maruyama step
        u = u + mu*dt + sqrt(2*sigma2*dt)*randn;

        % Threshold and reset
        if u >= theta
            if n > Ndiscard
                spikes = spikes + 1;
            end
            u = 0;
        end

        % Collect statistics after transient
        if n > Ndiscard
            count = count + 1;

            if u < 0
                below = below + 1;
            end

            if mod(n, sampleEvery) == 0
                samples(end+1) = u;
            end
        end
    end

    % Measured values
    Tmeas = T - Tdiscard;
    rate_meas = spikes / Tmeas;
    frac_meas = below / count;

    % Theoretical predictions
    rate_theory = mu/theta;
    frac_theory = sigma2/(mu*theta) * ...
        (1 - exp(-mu*theta/sigma2));

    fprintf('sigma^2 = %.2f\n', sigma2);
    fprintf('  firing rate: measured = %.4f, theory = %.4f\n', ...
        rate_meas, rate_theory);
    fprintf('  fraction u<0: measured = %.4f, theory = %.4f\n\n', ...
        frac_meas, frac_theory);

    % Histogram
    subplot(3,1,k);
    histogram(samples, 80, 'Normalization', 'pdf');
    hold on;

    % Theoretical stationary density
    x = linspace(min(samples), theta, 500);
    p = zeros(size(x));

    neg = x < 0;
    pos = x >= 0;

    p(neg) = (rate_theory/mu) * ...
        (1 - exp(-mu*theta/sigma2)) .* ...
        exp(mu*x(neg)/sigma2);

    p(pos) = (rate_theory/mu) * ...
        (1 - exp(mu*(x(pos)-theta)/sigma2));

    plot(x, p, 'LineWidth', 2);

    xlabel('u');
    ylabel('density');
    title(sprintf('\\sigma^2 = %.2f', sigma2));
    legend('simulation', 'theory', 'Location', 'northwest');
end