%% Explicit eulers
close all; clear all;

lw = 3;


%% Set solver parameters and simulation scenario
% ========================================================================
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
lambda = -1;
x0 = 1;
tspan = [0 10];
h = 0.1;

%% Our RK3 method

%% Compute numerical solution, analytical solution, and errors
% ========================================================================
[Tout,Xout,Eout] = RK3(@func,tspan,x0,h,lambda);
X = x0*exp(lambda*Tout);  % Real value
Eglobal = X-Xout;
EglobalK = Eglobal; % Normalized by the signal

Elocal = zeros(size(Eout));
ElocalK = zeros(size(Eout));  % Constant of error

%% Compute the actual local error

for i=2:length(Tout)
    Xlocal = Xout(i-1)*exp(lambda*(Tout(i)-Tout(i-1)));
    Elocal(i) = Xout(i) - Xlocal;
    ElocalK(i) = Elocal(i) / Xout(i-1);
end
% The first local error will be 0, becaute it is the initial conditions, so
% we get
Elocal(1) = Elocal(2);  
Eout(1) = Eout(2);  
%% Plotting
figure()

subplot(1,3,1); 
title('Real and Simulated Functions', 'fontSize',20,'fontWeight','Bold');
hold on
    plot(Tout, Xout, 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout, x0*exp(lambda*Tout), 'color',rand(1,3),'LineWidth',lw)  %e^(lambda*x)*x0
    legend('Simulated', 'Real')
    
    ylabel('U','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
    
hold off

subplot(1,3,2); title('Actual and Estimated Local Error');
hold on
    plot(Tout, log(abs(Eout)), 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout, log(abs(Elocal)), 'color',rand(1,3), 'LineWidth',lw)
    plot(Tout, log(abs(ElocalK)), 'color',rand(1,3), 'LineWidth',lw)
    
    legend('Estimated Local Error', 'Actual Local Error', 'Normalized Local Error')
    ylabel('log(|Error|)','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
    
hold off

subplot(1,3,3); title('Global error');
hold on
    plot(Tout, abs(x0*exp(lambda*Tout) - Xout), 'color',rand(1,3), 'LineWidth',lw)
    legend('Global Error')
    ylabel('|Error|','FontSize',12,'FontWeight','bold')
    xlabel('time','FontSize',12,'FontWeight','bold')
    
hold off


%% Study with the stepsize

N = 30;
tspan = [0 10];
hs = (1:N)/100;

LEs = zeros(N,1);
LTEs = zeros(N,1);
for ih = 1:N
    h = hs(1,ih);
    [Tout,Xout,Eout] = RK3(@func,tspan,x0,h,lambda);
    X = x0*exp(lambda*Tout);  % Real value
    Eglobal = X-Xout;

    Elocal = zeros(size(Eout));
    for i=2:length(Tout)
        Xlocal = Xout(i-1)*exp(lambda*(Tout(i)-Tout(i-1)));
        Elocal(i) = Xout(i) - Xlocal;
    end
    LEs(ih) = Elocal(2);
    LTEs(ih) = Eout(2);
end

%% Compute the actual local error
figure()

hsLog = log(hs);
LELog = log(abs(LEs));
LTELog = log(abs(LTEs));

slopeLE = (LELog(4) - LELog(3))/(hsLog(4) - hsLog(3));
slopeLTE = (LTELog(4) - LTELog(3))/(hsLog(4) - hsLog(3));
hold on
    title('Relation between local error and h','FontSize',12,'FontWeight','bold')
    plot(hsLog, LELog, 'color',rand(1,3), 'LineWidth',lw)
    plot(hsLog, LTELog, 'color',rand(1,3), 'LineWidth',lw)
    legend(strcat('Real. Slope = ',num2str(slopeLE)), strcat('Estimation. Slope = ',num2str(slopeLTE)))
    ylabel('log(|Error|)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
hold off


h = 0.1;
lambda = -1;

ts = 0:h:10;
N = size(ts);
N = N(2);
% Global error
Real = zeros(N,1);
Est = zeros(N,1);

for i = 1:N
    Real(i) = exp(lambda*0)*exp(lambda*ts(i));
    Est(i) = exp(lambda*0)*power (1 - h, i-1);
end

figure()
hold on
    plot(ts, abs(Real - Est))
    plot(ts, abs( Est))
    plot(ts, abs(Real))
    legend('Error', 'Estimation', 'Real')
 %   plot(ts, abs(Real - Est)./Est)
hold off
