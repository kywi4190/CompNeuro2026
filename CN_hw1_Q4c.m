clear; close all; clc

R = 1;
tm = 10;
ts = 2;
urest = 0;
uth = 1;

r = ts/tm;

tpeak = tm*ts/(ts-tm)*log(ts/tm);
Icrit = (uth-urest)/R * r^(1/(r-1));

I0s = [0.8 1 1.2]*Icrit;
t = linspace(0,40,4001);

peaks = zeros(1,3);

figure
hold on
for k = 1:3
    I0 = I0s(k);
    f = @(t,u) (-(u-urest) + R*I0*exp(-t/ts))/tm;
    [T,U] = ode45(f,t,urest);
    peaks(k) = max(U);
    plot(T,U)
end
yline(uth,'--')
xlabel('t (ms)')
ylabel('u (mV)')
legend('Below','At','Above','Threshold')
hold off

fprintf('t_peak = %.4f ms\n',tpeak)
fprintf('I_0^crit = %.4f mA\n\n',Icrit)

fprintf('I_0 (mA)       Peak u (mV)\n')
for k = 1:3
    fprintf('%.4f       %.4f\n',I0s(k),peaks(k))
end

% Icrit versus tau_s
ts_vals = logspace(log10(0.5),log10(100),500);
r = ts_vals/tm;

x = log(r)./(r-1);
x(abs(r-1) < 1e-8) = 1;

Icrit_vals = (uth-urest)/R .* exp(x);

figure
loglog(ts_vals,Icrit_vals)
xlabel('\tau_s (ms)')
ylabel('I_0^{crit} (mA)')
grid on

% Written answer:
% (c)
% as ts -> 0, Icrit -> inf
% as ts -> inf, Icrit -> (uth - urest)/R