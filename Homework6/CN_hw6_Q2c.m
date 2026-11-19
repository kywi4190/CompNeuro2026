clear; clc; close all;

N = 4e5;
eta = 0.02;
xs = 8/17;
ys = 9/17;
P = [0.40 0.60; 0.40 0.25];

% (i) Mean flow from part (b).
flow = @(t,z) 0.85*eta * [(z(1)*(1-z(1)))^2*(z(2)-ys); -(z(2)*(1-z(2)))^2*(z(1)-xs)];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,z] = ode45(flow,0:10:40000,[0.2; 0.3],opts);

K = zeros(size(t));
for k = 1:length(t)
    K(k) = integral(@(s) (s-xs)./(s.*(1-s)).^2,xs,z(k,1)) + integral(@(s) (s-ys)./(s.*(1-s)).^2,ys,z(k,2));
end

fprintf('Mean-flow x range: [%.5f, %.5f]\n',min(z(:,1)),max(z(:,1)));
fprintf('Mean-flow y range: [%.5f, %.5f]\n',min(z(:,2)),max(z(:,2)));
fprintf('Variation in K (max-min): %.3e\n',max(K)-min(K));

figure;
subplot(1,3,1);
plot(z(:,1),z(:,2),'b'); hold on;
plot(xs,ys,'rx','LineWidth',2);
xlabel('x'); ylabel('y'); title('Mean flow');
axis([0 1 0 1]); axis square; grid on;

% (ii) Random games with and without the entropy bonus
taus = [0 0.05];
early = 1:N/4;
late = 3*N/4+1:N;

for j = 1:length(taus)
    rng(1);
    tau = taus(j);
    theta = [0.3; -0.3];
    b = [0; 0];
    XY = zeros(N,2);

    for n = 1:N
        p = 1./(1+exp(-theta));
        XY(n,:) = p';

        % Action 1 means shoot 3 or guard the paint.
        a = rand(2,1) < p;
        made = rand < P(1+a(1),2-a(2));
        points = (2+a(1))*made;
        reward = [points; -points];

        % Update both players then the running-average baseline
        theta = theta + eta*((reward-b).*(a-p) - tau*theta.*p.*(1-p));
        b = b + (reward-b)/n;
    end

    avg = mean(XY(late,:),1);
    d = sqrt((XY(:,1)-xs).^2 + (XY(:,2)-ys).^2);
    fprintf('\ntau = %.2f\n',tau);
    fprintf('Last-quarter mean (x,y): (%.5f, %.5f)\n',avg(1),avg(2));
    fprintf('Mean distance, first quarter: %.5f\n',mean(d(early)));
    fprintf('Mean distance, last quarter:  %.5f\n',mean(d(late)));

    subplot(1,3,j+1);
    plot(XY(:,1),XY(:,2),'b'); hold on;
    plot(xs,ys,'rx','LineWidth',2);
    xlabel('x'); ylabel('y'); title(sprintf('\\tau = %.2f',tau));
    axis([0 1 0 1]); axis square; grid on;
end