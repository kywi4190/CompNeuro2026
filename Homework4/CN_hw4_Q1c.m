clear; close all; clc;
rng(1);

N = 200;
sigma_e = 1;
sigma_i_values = [2, 3];
dt = 0.02;
% Time unit: tau_w / (q_S - q_D)
T = 40;
nTrials = 20;

a = (0:N-1)';
% Periodic distance from site 0
d = min(a, N-a);
periods = zeros(nTrials, 2);

for j = 1:2
    sigma_i = sigma_i_values(j);

    K = exp(-d.^2/(2*sigma_e^2)) ...
        - (sigma_e/sigma_i)*exp(-d.^2/(2*sigma_i^2));
    Khat = fft(K);

    for trial = 1:nTrials
        % w represents w_-
        w = 0.01*randn(N, 1);

        for step = 1:round(T/dt)
            % Forward Euler, using circular convolution
            w = w + dt*real(ifft(Khat .* fft(w)));
            w = max(-1, min(1, w));
        end

        % Two sign changes per stripe pair, including the ring boundary
        right = (w > 0);
        nPairs = sum(right ~= right([2:N, 1]))/2;
        periods(trial, j) = N/nPairs;

        if trial == 1
            finalPattern = w;
        end
    end

    % Continuum prediction
    kstar = sqrt(4*log(sigma_i/sigma_e)/(sigma_i^2-sigma_e^2));
    continuumPeriod = 2*pi/kstar;

    % Allowed ring modes: k_m = 2*pi*m/N, m = 1,...,N/2.
    % Their growth rates are the corresponding FFT entries
    [~, mstar] = max(real(Khat(2:N/2+1)));
    ringPeriod = N/mstar;

    fprintf('\nsigma_i = %g; periods in site spacings\n', sigma_i);
    fprintf('Measured: %.4g +/- %.4g (mean +/- SD, 20 runs)\n', ...
        mean(periods(:,j)), std(periods(:,j)));
    fprintf('Continuum: %.4g\n', continuumPeriod);
    fprintf('Fastest ring mode: %.4g (m = %d)\n', ringPeriod, mstar);

    % Plot one final pattern, closing the ring at a = N
    figure;
    plot([a; N], [finalPattern; finalPattern(1)], 'LineWidth', 1.2);
    xlabel('Position a (site spacings)');
    ylabel('Eye preference w_-');
    title(sprintf('N = 200, \\sigma_e = 1, \\sigma_i = %g, t = 40', ...
        sigma_i));
    xlim([0, N]);
    ylim([-1.1, 1.1]);
    grid on;
end
