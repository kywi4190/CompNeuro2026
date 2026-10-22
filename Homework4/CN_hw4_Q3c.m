clear; close all; clc;
rng(1);

N = 500;
alpha = 0.02:0.04:0.62;
numCues = 10;
unstableFraction = zeros(size(alpha));
meanOverlap = zeros(size(alpha));

for k = 1:numel(alpha)
    P = round(alpha(k)*N);
    X = 2*randi([0 1],N,P) - 1;

    J = X*X';
    % J = N*W: same signs, exact zero fields
    J(1:N+1:end) = 0;

    % Fraction of unstable bits over all stored patterns
    unstable = X.*(J*X) < 0;
    unstableFraction(k) = mean(unstable(:));

    cueIndices = randperm(P,numCues);
    overlap = zeros(numCues,1);

    for trial = 1:numCues
        target = X(:,cueIndices(trial));
        s = target;

        % Flip exactly 10% of the bits
        flipped = randperm(N,round(0.10*N));
        s(flipped) = -s(flipped);

        % Update one unit at a time until a whole sweep changes nothing
        changed = true;
        while changed
            changed = false;
            for i = randperm(N)
                % sgn(0) = +1
                newSpin = 2*(J(i,:)*s >= 0) - 1;
                if newSpin ~= s(i)
                    s(i) = newSpin;
                    changed = true;
                end
            end
        end

        overlap(trial) = target'*s/N;
    end

    meanOverlap(k) = mean(overlap);
end

theory = 0.5*erfc(1./sqrt(2*alpha));

figure('Color','w','Position',[100 100 1000 400]);

subplot(1,2,1);
plot(alpha,unstableFraction,'bo-',alpha,theory,'k--','LineWidth',1.5);
xlabel('\alpha = P/N');
ylabel('Fraction of unstable bits');
title('Stored-pattern stability: N = 500');
legend('Simulation','Gaussian prediction','Location','northwest');
grid on;

subplot(1,2,2);
plot(alpha,meanOverlap,'bo-','LineWidth',1.5);
xlabel('\alpha = P/N');
ylabel('Mean overlap m');
title('Recall: 10 cues, 10% of bits flipped');
grid on;

saveas(gcf,'q3c_hopfield.png');