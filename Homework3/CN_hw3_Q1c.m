clear; clc; close all;

%% Parameters
eps = 0.5;
Ivals = linspace(0.9, 3, 700);

% transient time
Ttrans = 200;
% time over which spikes are counted
Tmeasure = 400;
% step used to bracket threshold crossings
dt_search = 0.01;

% From part (a)
Istar = 1/(1-exp(-1));
tongueWidth = eps/sqrt(1 + 4*pi^2);

Ileft  = Istar - tongueWidth;
Iright = Istar + tongueWidth;

%% Simulate forced system
firingNum = zeros(size(Ivals));

for j = 1:length(Ivals)

    I = Ivals(j);

    % Treat t = 0 as a reset: v(0) = 0
    Tn = 0;
    spikeCount = 0;

    while Tn < Ttrans + Tmeasure

        % Find next spike
        Tnext = nextSpike(Tn, I, eps, dt_search);

        if isinf(Tnext) || Tnext > Ttrans + Tmeasure
            break;
        end

        % Count only after transient
        if Tnext >= Ttrans
            spikeCount = spikeCount + 1;
        end

        Tn = Tnext;
    end

    % forcing period = 1, so this is spikes per period
    firingNum(j) = spikeCount/Tmeasure;
end


%% Unforced firing number
unforced = zeros(size(Ivals));

idx = Ivals > 1;
unforced(idx) = 1 ./ log(Ivals(idx)./(Ivals(idx)-1));


%% Plot
figure;
plot(Ivals, firingNum, 'b', 'LineWidth', 1.5);
hold on;
plot(Ivals, unforced, 'k--', 'LineWidth', 1.5);

% 1:1 tongue boundaries
xline(Ileft,  'r--', 'LineWidth', 1.2);
xline(Iright, 'r--', 'LineWidth', 1.2);

% Show firing number = 1
yline(1, ':');

xlabel('I');
ylabel('Firing number (spikes per forcing period)');
title('\epsilon = 0.5');
legend('Forced', 'Unforced', ...
       '1:1 tongue boundary', '1:1 tongue boundary', ...
       'Location', 'southeast');
grid on;

fprintf('Predicted 1:1 tongue: %.4f <= I <= %.4f\n', ...
        Ileft, Iright);


%% Find spike phase at I = 1.6
I = 1.6;

Tn = 0;
phases = [];

while Tn < Ttrans + 20

    Tnext = nextSpike(Tn, I, eps, dt_search);

    if isinf(Tnext)
        break;
    end

    if Tnext > Ttrans
        % Phase within forcing period
        phases(end+1) = mod(Tnext,1);
    end

    Tn = Tnext;
end

fprintf('Spike phase at I = 1.6: phi = %.5f\n', mean(phases));


% ---------------------------------------------------
%  Function for finding the next threshold crossing
function Tnext = nextSpike(Tn, I, eps, dt)

    G = @(t) I + eps/(1 + 4*pi^2) .* ...
        (sin(2*pi*t) - 2*pi*cos(2*pi*t));

    % Exact voltage after reset at Tn:
    % v(t) = G(t) - G(Tn)*exp(-(t-Tn))
    v = @(t) G(t) - G(Tn).*exp(-(t-Tn));

    % We want the first root of v(t) - 1
    f = @(t) v(t) - 1;

    t1 = Tn;
    f1 = f(t1);

    % Search forward until threshold is crossed
    maxWait = 30;

    for t2 = Tn+dt : dt : Tn+maxWait

        f2 = f(t2);

        if f1 < 0 && f2 >= 0
            Tnext = fzero(f, [t1 t2]);
            return;
        end

        t1 = t2;
        f1 = f2;
    end

    % No spike found
    Tnext = inf;
end
