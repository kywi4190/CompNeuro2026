
clear; clc; close all;

L = 200;
phi0 = 0.4;
epsVals = [0.005, 0.02, 0.05];

for e = epsVals
    tEnd = 5/e;
    t = 0;

    % Neuron 1 starts at reset.
    % Neuron 2 starts 0.4 cycles behind.
    v = [-L; -cot(pi*(1-phi0))];

    % Columns: [spike time, neuron number]
    spikes = [];

    % Coupled QIF equations
    f = @(t,v) [1 + v(1)^2 + e*(v(2)-v(1));
        1 + v(2)^2 + e*(v(1)-v(2))];

    opts = odeset('Events', @(t,v) spikeEvent(t,v,L), ...
        'RelTol', 1e-8, 'AbsTol', 1e-9, ...
        'MaxStep', 0.1);

    % Simulate, resetting neurons whenever they spike
    while t < tEnd
        [~,~,te,ve,ie] = ode45(f, [t tEnd], v, opts);

        if isempty(te)
            break;
        end

        t = te(end);
        v = ve(end,:)';
        v(ie) = -L;

        spikes = [spikes; te(:), ie(:)];
    end

    % Extract spike times
    t1 = spikes(spikes(:,2)==1, 1);
    t2 = spikes(spikes(:,2)==2, 1);

    % Measure phase lag once per cycle
    phi = [];
    times = [];

    for k = 1:length(t1)-1
        j = find(t2 > t1(k), 1);

        if isempty(j)
            continue;
        end

        lag = (t2(j)-t1(k))/(t1(k+1)-t1(k));

        if lag > 0 && lag < 0.5
            phi(end+1) = lag;
            times(end+1) = t1(k);
        end
    end

    % Fit log(tan(pi*phi)) = constant - rate*t
    % Exclude early transients and very small lags
    idx = phi > 0.02 & phi < 0.35;

    p = polyfit(times(idx), log(tan(pi*phi(idx))), 1);
    rate = -p(1);

    fprintf('epsilon = %.3f | fitted = %.5f | predicted = %.5f\n', ...
        e, rate, 2*e);

    % Compare measured lag with analytical prediction
    if e == 0.02
        tt = linspace(0, tEnd, 500);
        prediction = atan(tan(pi*phi0)*exp(-2*e*tt))/pi;

        figure;
        plot(times, phi, 'w.', 'MarkerSize', 10);
        hold on;
        plot(tt, prediction, 'r-', 'LineWidth', 1.5);

        xlabel('Time');
        ylabel('Phase lag \phi');
        title('Gap-junction QIF neurons (\epsilon = 0.02)');
        legend('Measured lag', 'Phase model');
        grid on;
    end
end

% Stop integration when either neuron reaches +L
function [value, isterminal, direction] = spikeEvent(~,v,L)
value = v - L;
isterminal = [1; 1];
direction = [1; 1];
end
