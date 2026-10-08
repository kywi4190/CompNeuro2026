clear; close all; clc;

% size of voltage kick
delta = 1e-3;
% phase grid
phi = linspace(0,0.99,100);


% ============================================================
% 1. LIF MODEL

tau = 1;
vr  = 0;
vth = 1;
I   = 1.5;

% Period from part (a)
Delta = tau*log((I*tau - vr)/(I*tau - vth));

R_num = zeros(size(phi));

for k = 1:length(phi)

    % Time and voltage at phase phi
    tk = phi(k)*Delta;
    vk = I*tau + (vr - I*tau)*exp(-tk/tau);

    % Apply voltage kick
    vk = vk + delta;

    % Find time from kick until threshold
    opts = odeset('Events', @(t,v) threshold(t,v,vth), ...
                  'RelTol',1e-9,'AbsTol',1e-11);

    [~,~,te] = ode45(@(t,v) I - v/tau, [0 2*Delta], vk, opts);

    % Total time from reset to perturbed spike
    Tprime = tk + te(1);

    % Numerically measured PRC
    R_num(k) = (Delta - Tprime)/(Delta*delta);
end

% Exact formula from part (a)
R_exact = tau*exp(phi*Delta/tau) ./ (Delta*(I*tau - vr));

figure;
plot(phi,R_num,'o','MarkerSize',4);
hold on;
plot(phi,R_exact,'LineWidth',2);
xlabel('\phi');
ylabel('R(\phi)');
title('LIF Phase Response Curve');
legend('Direct perturbation','Analytic PRC','Location','best');
grid on;

fprintf('LIF\n');
fprintf('Period = %.5f\n',Delta);
fprintf('R minimum = %.5f\n',min(R_exact));
fprintf('R maximum = %.5f\n',max(R_exact));
fprintf('R is never negative -> Type I PRC\n\n');


% ============================================================
% 2. FITZHUGH-NAGUMO MODEL

FHN = @(t,y) [ ...
    y(1) - y(1)^3/3 - y(2) + 0.5;
    0.08*(y(1) + 0.7 - 0.8*y(2)) ];

% First run for a long time so transients die away
y0 = [-1; -0.5];

optsCross = odeset('Events',@upwardCross, ...
                   'RelTol',1e-9,'AbsTol',1e-11, ...
                   'MaxStep',0.1);

[~,~,tCross,yCross] = ode45(FHN,[0 500],y0,optsCross);

% Period from final two upward v = 0 crossings
DeltaF = tCross(end) - tCross(end-1);

% Start exactly at a late, stable upward crossing
yStart = yCross(end,:)';

% Find the limit-cycle state at each desired phase
tPhase = phi*DeltaF;

opts = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.1);

[~,Yphase] = ode45(FHN,tPhase,yStart,opts);

RF = zeros(size(phi));

% Measure crossing several cycles later
nCross = 4;

for k = 1:length(phi)

    % State on unperturbed limit cycle
    yKick = Yphase(k,:)';

    % Kick only v
    yKick(1) = yKick(1) + delta;

    % Integrate perturbed trajectory and detect upward crossings
    [~,~,te] = ode45(FHN,[0 (nCross+1)*DeltaF], ...
                       yKick,optsCross);

    % Absolute time of nth future crossing
    Tperturbed = phi(k)*DeltaF + te(nCross);

    % Without perturbation, that crossing would occur here
    Tunperturbed = nCross*DeltaF;

    % PRC
    RF(k) = (Tunperturbed - Tperturbed)/(DeltaF*delta);
end

figure;
plot(phi,RF,'o-','MarkerSize',4);
yline(0,'k--');
xlabel('\phi');
ylabel('R(\phi)');
title('FitzHugh-Nagumo Phase Response Curve');
grid on;

% Requested quantities
[Rmin,imin] = min(RF);
[Rmax,imax] = max(RF);

neg = phi(RF < 0);

fprintf('FitzHugh-Nagumo\n');
fprintf('Period = %.5f\n',DeltaF);
fprintf('R minimum = %.5f at phi = %.3f\n',Rmin,phi(imin));
fprintf('R maximum = %.5f at phi = %.3f\n',Rmax,phi(imax));

if ~isempty(neg)
    fprintf('R < 0 approximately for %.3f < phi < %.3f\n', ...
            min(neg),max(neg));
end


% ============================================================
% Event functions

function [value,isterminal,direction] = threshold(~,v,vth)
    value = v - vth;
    isterminal = 1;
    direction = 1;
end

function [value,isterminal,direction] = upwardCross(~,y)
    % v = 0
    value = y(1);
    % keep integrating
    isterminal = 0;
    % upward crossings only
    direction = 1;
end