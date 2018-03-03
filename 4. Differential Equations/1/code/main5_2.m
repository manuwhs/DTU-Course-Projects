%% Explicit eulers
close all; clear all;

lw = 3;

limits = 10;
alpha = -limits:0.1:limits;
beta = -limits:0.1:limits;

gamma = 1-1/sqrt(2);
a31 = (1-gamma)/2;
AT = [0 gamma a31;0 gamma a31;0 0 gamma];
c  = [0; 2*gamma; 1];
b  = AT(:,3);
bhat = [    (6*gamma-1)/(12*gamma); ...
    1/(12*gamma*(1-2*gamma)); ...
    (1-3*gamma)/(3*(1-2*gamma))    ];
d  = b-bhat;
p  = 2;
phat = 3;
s = 3;

A = AT.';
AT = A.';

nreal = length(alpha);
nimag = length(beta);
I = eye(size(A));
e = ones(size(A,1),1);

for kreal = 1:nreal
    for kimag = 1:nimag
        z = alpha(kreal) + 1i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b.'*tmp;
        Ehat = z*d.'*tmp;
        f = exp(z); % This is the real variation
        E = R-f;    % This is the estimated variation
        EhatmE = Ehat-E;
        
        absR(kimag,kreal) = abs(R);
        absEhatmE(kimag,kreal) = abs(EhatmE);
        absEhat(kimag,kreal) = abs(Ehat);
        absE(kimag,kreal) = abs(E);
        absF(kimag,kreal) = abs(f);
    end
end


figure(2)
fs = 14;
subplot(221)
imagesc(alpha,beta,absR,[0 1]); % In the plotting we show up until maximum
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',fs);
ylabel('imag','fontsize',fs);
title('|R(z)|','fontsize',fs)

subplot(222)
imagesc(alpha,beta,absE,[0 1]); % In the plotting we show up until maximum
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',fs);
ylabel('imag','fontsize',fs);
title('|E(z)|','fontsize',fs)

subplot(223)
imagesc(alpha,beta,absEhat,[0 1]); % In the plotting we show up until maximum
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',fs);
ylabel('imag','fontsize',fs);
title('|Ehat(z)|','fontsize',fs)

subplot(224)
imagesc(alpha,beta,absEhatmE,[0 1]); % In the plotting we show up until maximum
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',fs);
ylabel('imag','fontsize',fs);
title('|Ehat(z) - E(z)|','fontsize',fs)

