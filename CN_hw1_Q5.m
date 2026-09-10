clear; close all; clc

%% (a) Rheobase, gK = 36
rh1 = rheobase(36);

fprintf('Rheobase (gK = 36): %.3f uA/cm^2\n',rh1)

[t1,V1] = hh_sim(rh1-0.01);
[t2,V2] = hh_sim(rh1+0.01);

yl = [min([V1;V2])-5 max([V1;V2])+5];

figure
subplot(1,2,1)
plot(t1,V1)
xlim([0 300]); ylim(yl)
xlabel('Time (ms)')
ylabel('V (mV)')
title('Just below rheobase')

subplot(1,2,2)
plot(t2,V2)
xlim([0 300]); ylim(yl)
xlabel('Time (ms)')
ylabel('V (mV)')
title('Just above rheobase')

%% (b) f-I curve
I1 = linspace(rh1+0.01,2*rh1,25);
f1 = zeros(size(I1));

for k = 1:length(I1)
    [~,~,~,f1(k)] = hh_sim(I1(k));
end

fprintf('Onset rate: %.3f Hz\n',f1(1))

figure
plot(I1,f1,'o-')
xlabel('I (\muA/cm^2)')
ylabel('Firing rate (Hz)')
grid on

%% (c) Change gK to 30
rh2 = rheobase(30);

fprintf('\n gK      Rheobase\n')
fprintf(' 36        %.3f\n',rh1)
fprintf(' 30        %.3f\n',rh2)

I2 = linspace(rh2+0.01,2*rh1,25);
f2 = zeros(size(I2));

for k = 1:length(I2)
    [~,~,~,f2(k)] = hh_sim(I2(k),[],30);
end

figure
plot(I1,f1,'o-')
hold on
plot(I2,f2,'o-')
xline(rh1,'--')
xline(rh2,'--')
xlabel('I (\muA/cm^2)')
ylabel('Firing rate (Hz)')
legend('g_K = 36','g_K = 30','Rheobase 36','Rheobase 30')
grid on
hold off

%% Bisection
function rh = rheobase(gK)

lo = 0;
hi = 10;

while hi-lo > 1e-4
    mid = (lo+hi)/2;
    [~,~,spk] = hh_sim(mid,[],gK);

    if length(spk) >= 2
        hi = mid;
    else
        lo = mid;
    end
end

rh = (lo+hi)/2;
end

% Written answers:
% (a)
% Repetitive firing criterion is at least 2 spikes after the startup.

% (b)
% The firing rate switches on at a nonzero value, so Hodgkin-Huxley
% resembles the Hopf bifurcation rather than the saddle-node bifurcation.
% Physiologically, this means a non-zero firing rate when repetitive firing begins.

% (c)
% Lowering gK lowered the rheobase because there is less outward K+
% current opposing depolarization. The curve shifted left and its
% shape appeared to change somewhat.