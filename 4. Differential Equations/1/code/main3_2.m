%% Explicit eulers
close all; clear all;

lw = 3;

alpha = -5:0.01:5;
beta = -5:0.01:5;

A = [0 0 0;
    1/2 0 0;
    -1 2 0];
b = [1/6; 2/3; 1/6];
c = [0; 1; 1/2];
d = [1/12; -1/6; 1/12];
s = 3;

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
        
        Ehat = z*d.'*tmp; % This is the estimated local error
        f = exp(z); % This is the real variation
        E = R-f;    % This is the real local error
        EhatmE = Ehat-E; % This is the difference between the estimated error and the real error.
        
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

