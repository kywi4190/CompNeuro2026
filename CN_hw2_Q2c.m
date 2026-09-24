clear; clc; close all;
V = linspace(-90,0,1000);
V0 = 16.13;

hold on
for Mg = [0 0.1 1]
    B = 1 ./ (1 + (Mg/3.57).*exp(-V/V0));
    % EN = 0 mV
    % g_bar not specified, set it to 1
    IN = B .* V;
    plot(V,IN,'DisplayName',sprintf('[Mg^{2+}] = %.1f mM',Mg))
end

xlabel('V (mV)')
ylabel('I_N / g_bar','Interpreter','tex')
title('NMDA Current vs. Voltage')
legend('Location','southeast')
grid on