%% Explicit eulers
close all; clear all;

tspan = 10;
Nsteps = 20;
lw = 3;
[tnList_0,ynList_0] = ExplicitEulers(@func,tspan,Nsteps*2,1);
[tnList_1,ynList_1] = ExplicitTrapezoid(@func,tspan,Nsteps,1);
[tnList_2,ynList_2] = Midpoint(@func,tspan,Nsteps,1);
[tnList_4,ynList_4] = RK4(@func,tspan,Nsteps,1);

tspan = 50;
Y0 = [2 ; 2];


%% Set solver parameters and simulation scenario
% ========================================================================
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
lambda = -1;
x0 = 1;
tspan = [0 10];
h = 0.25;

%% Compute numerical solution, analytical solution, and errors
% ========================================================================
[Tout,Xout,Eout] = ...
ERKSolverErrorEstimation(@func,tspan,x0,h,solver,lambda);
X = x0*exp(lambda*Tout);
Eglobal = X-Xout;
Elocal = zeros(size(Eout));

%% Compute the bounded local error
for i=2:length(Tout)
    Xlocal = Xout(i-1)*exp(lambda*(Tout(i)-Tout(i-1)));
    Elocal(i) = Xout(i) - Xlocal;
end

figure()
subplot(2,2,1); title('Real and Estimation');
hold on
    plot(Tout, Xout, 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout, exp(1*Tout)*1, 'color',rand(1,3),'LineWidth',lw)  %e^(lambda*x)*x0
    legend('DOPRI', 'real')
hold off

subplot(2,2,2); title('Local Error');
hold on
    plot(Tout, Eout, 'color',rand(1,3), 'LineWidth',lw)
    legend('Actual error')
hold off

subplot(2,2,3); title('Bounded Error');
hold on
    plot(Tout, Elocal, 'color',rand(1,3), 'LineWidth',lw)
    legend('Bounded error')
hold off

%% Impicit Eulers