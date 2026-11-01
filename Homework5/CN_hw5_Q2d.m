clear; clc; close all

beta = 1;
g = 2;
tau_a = 120;
eta = 1000;
IR = 1.3;

% Each row gives w and I_L.
cases = [0.8 1.3; 0.8 1.45; 0.8 1.6; 0.3 1.3];
f = @(x) 1./(1 + exp(-eta*x));
options = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
    'MaxStep',1,'Events',@handover);

fprintf('(ii) w = 0.8, I_R = 1.3\n');
fprintf('I_L    TL sim  TL theory   TR sim  TR theory   TL error(%%)  TR error(%%)\n');

for k = 1:size(cases,1)
    w = cases(k,1);
    IL = cases(k,2);

    rhs = @(t,y) [
        -y(1) + f(w*y(1) - beta*y(2) - g*y(3) + IL);
        -y(2) + f(w*y(2) - beta*y(1) - g*y(4) + IR);
        (y(1) - y(3))/tau_a;
        (y(2) - y(4))/tau_a];

    [t,y,te,~,ie] = ode15s(rhs,[0 30*tau_a],[1;0;0;0],options);

    % Count only switches to the opposite percept.
    keep = diff([1; ie]) ~= 0;
    te = te(keep);
    ie = ie(keep);
    T = diff(te);
    winner = ie(1:end-1);

    % Remove the first three complete dominance intervals.
    T = T(4:end);
    winner = winner(4:end);
    TL = mean(T(winner == 1));
    TR = mean(T(winner == 2));

    if k <= 3
        BL = (IL-beta)/g;
        BR = (IR-beta)/g;
        TL_theory = tau_a*log((1-BL)/BR);
        TR_theory = tau_a*log((1-BR)/BL);
        errL = 100*(TL-TL_theory)/TL_theory;
        errR = 100*(TR-TR_theory)/TR_theory;

        fprintf('%.2f %8.3g %10.3g %8.3g %10.3g %+12.3g %+12.3g\n', ...
            IL,TL,TL_theory,TR,TR_theory,errL,errR);
    else
        Tdom = mean(T);
        T_theory = tau_a*log((w+IL)/(g-w-IL));

        fprintf('\n(iii) w = 0.3, I_L = I_R = 1.3\n');
        fprintf('Tdom = %.3g, theory = %.3g\n',Tdom,T_theory);
        fprintf('Difference = %+.3g time units (%+.3g%%)\n', ...
            Tdom-T_theory,100*(Tdom-T_theory)/T_theory);
    end

    if k == 1
        B = (IR-beta)/g;
        figure;

        subplot(2,1,1)
        plot(t,y(:,1:2),'LineWidth',1);
        ylabel('Rate');
        legend('u_L','u_R','Location','best');
        grid on

        subplot(2,1,2)
        p = plot(t,y(:,3:4),'LineWidth',1);
        hold on
        yline(B,'k--','B');
        yline(1-B,'k:','1-B');
        xlabel('Time (rate time constants)');
        ylabel('Adaptation');
        legend(p,{'a_L','a_R'},'Location','best');
        grid on

        sgtitle({'I_L=I_R=1.3, w=0.8, \beta=1, g=2', ...
            '\tau_a=120, \eta=1000'});
    end
end

% Each event marks a percept becoming dominant.
function [value,isterminal,direction] = handover(~,y)
value = [y(1)-y(2)-0.5; y(2)-y(1)-0.5];
isterminal = [0;0];
direction = [1;1];
end