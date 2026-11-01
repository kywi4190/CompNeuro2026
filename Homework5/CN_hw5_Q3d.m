clear; close all; clc;

dx = 0.02;
dt = 0.02;
x = -20:dx:20;
t = (0:dt:60)';
N = length(x);

W = @(z) 2*sign(z).*(exp(-abs(z)/2) - exp(-abs(z)));

% Pad the FFT for the full convolution.
lags = (1-N:N-1)*dx;
w = 2*exp(-abs(lags)) - exp(-abs(lags)/2);
M = 2^nextpow2(3*N - 2);
K = fft(w, M);

% The first four runs test stability; the last three measure rates.
kappas = [0.3 0.3 0.3 0.3 0.1 0.3 0.45];
offsets = [-0.02 0.02 -0.05 0.05 0.05 0.05 0.05];
a = zeros(length(t), 7);

for j = 1:7
    kappa = kappas(j);
    S = sqrt(1 - 2*kappa);
    aminus = -log((1 + S)/2);
    aplus = -log((1 - S)/2);

    a0 = aplus + offsets(j);
    if j <= 2
        a0 = aminus + offsets(j);
    end
    u = W(x + a0) - W(x - a0);

    for n = 1:length(t)
        % Interpolate the right edge; no active region means width zero.
        r = find(u >= kappa, 1, 'last');
        if ~isempty(r)
            a(n,j) = x(r) + dx*(kappa-u(r))/(u(r+1)-u(r));
        end

        % Take one Euler step using the padded convolution.
        f = 1./(1 + exp(-100*(u-kappa)));
        c = real(ifft(K.*fft(f, M)));
        u = u + dt*(-u + dx*c(N:2*N-1));
    end
end

figure;
h = plot(t, a(:,1:4), 'LineWidth', 1.4);
hold on;
S = sqrt(1 - 2*0.3);
yline(-log((1+S)/2), ':', 'a_-');
yline(-log((1-S)/2), ':', 'a_+');
xlabel('Time t');
ylabel('Active half-width a(t)');
title('\kappa = 0.3, \Delta x = \Delta t = 0.02, sigmoid slope = 100');
legend(h, {'a_0 = a_- - 0.02', 'a_0 = a_- + 0.02', ...
    'a_0 = a_+ - 0.05', 'a_0 = a_+ + 0.05'}, ...
    'Location', 'best');
grid on;

figure;
fprintf('kappa    fitted rate    theory\n');

for j = 5:7
    kappa = kappas(j);
    S = sqrt(1 - 2*kappa);
    aplus = -log((1 - S)/2);
    theory = -S*(1-S)/(1 + S*(1-S)/2);
    err = abs(a(:,j) - aplus);

    % Fit before the small final width error becomes important.
    fit = t >= 2 & t <= 10;
    p = polyfit(t(fit), log(err(fit)), 1);

    fprintf('%.2g      %.4g        %.4g\n', kappa, p(1), theory);

    subplot(1,3,j-4);
    semilogy(t, err, 'LineWidth', 1.4);
    hold on;
    semilogy(t(fit), exp(polyval(p,t(fit))), '--', 'LineWidth', 1.4);
    xlabel('Time t');
    ylabel('|a(t) - a_+|');
    title(sprintf('\\kappa = %.2g', kappa));
    legend('Simulation', 'Exponential fit', 'Location', 'north');
    grid on;
end
sgtitle('a_0 = a_+ + 0.05, \Delta x = \Delta t = 0.02, sigmoid slope = 100');