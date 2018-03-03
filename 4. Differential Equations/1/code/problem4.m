%% Question 1: mu = 3 P-controller
close all; clear all;

% Paramters
tspan = [0 15];
N = 1000;
y0 = [1;1];
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
h = (tspan(2)-tspan(1))/N;
abstol=10e-7;
reltol = 10e-7;

% Solvers
[tnList1,ynList1,hList1,rList1,nfun1] = ExplicitEulersAdaptiveStep(...
    @VanDerPol,tspan,N,y0,abstol,reltol,'P');
[tnList2,ynList2,hList2,rList2,nfun2] = ImplicitEulersAdaptiveStep(...
    @VanDerPol,@jacVanDerPol,tspan,N,y0,abstol,reltol,'P');
[tnList3,ynList3,hList3,rList3,nfun3] = RK4AdaptiveStep(...
    @VanDerPol,tspan,N,y0,abstol,reltol);
[tnList4,ynList4,hList4,rList4,nfun4] = ImplicitTrapezoidAdaptiveStep(...
    @VanDerPol,@jacVanDerPol,tspan,N,y0,abstol,reltol,'P');
% [tnList,ynList,Eout] = ERKSolverAdaptiveStep(...
%     @VanDerPol,tspan,y0,h,solver,abstol,reltol);

% Plotting
figure(1)
lw = 2;
subplot(2,1,1)
hold on
    plot(tnList1(1,:),ynList1(1,:),'-.','LineWidth',lw)
    plot(tnList2(1,:),ynList2(1,:),'-.','LineWidth',lw+1)
    plot(tnList3(1,:),ynList3(1,:),':','LineWidth',lw)
    plot(tnList4(1,:),ynList4(1,:),':','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4');
    set(leg,'FontSize',14);
    ylabel('x','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

subplot(2,1,2)
hold on
    plot(tnList1(1,:),ynList1(2,:),'-.','LineWidth',lw)
    plot(tnList2(1,:),ynList2(2,:),'-.','LineWidth',lw+1)
    plot(tnList3(1,:),ynList3(2,:),':','LineWidth',lw)
    plot(tnList4(1,:),ynList4(2,:),':','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4');
    set(leg,'FontSize',14);
    ylabel('y','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

figure(2)
lw = 2;
subplot(2,1,1)
hold on
    plot(tnList1(2:end),hList1,'-.','LineWidth',lw)
    plot(tnList2(2:end),hList2,'-.','LineWidth',lw)
    plot(tnList3(2:end),hList3,':','LineWidth',lw)
    plot(tnList4(2:end),hList4,'--','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4');
    set(leg,'FontSize',14);
    ylabel('h','FontSize',14,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',14,'FontWeight','bold')
    grid on
hold off

subplot(2,1,2)
ve = ones(1,length(tnList1)-1);
hold on
    plot(tnList1(2:end),rList1,'-.','LineWidth',lw)
    plot(tnList2(2:end),rList2,'-.','LineWidth',lw)
    plot(tnList3(2:end),rList3,':','LineWidth',lw)
    plot(tnList4(2:end),rList4,'--','LineWidth',lw)
    plot(tnList1(2:end),ve,'LineWidth',3)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4','MaxTolerance');
    set(leg,'FontSize',14);
    ylim([0 1.3])
    ylabel('r','FontSize',14,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',14,'FontWeight','bold')
    grid on
hold off

%% Question 2: mu = 3 PI
close all; clear all;

% Paramters
tspan = [0 15];
N = 1000;
y0 = [1;1];
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
h = (tspan(2)-tspan(1))/N;
abstol=10e-5;
reltol = 10e-5;

% Solvers
[tnList1,ynList1,hList1,rList1,nfun1] = ExplicitEulersAdaptiveStep(...
    @VanDerPol,tspan,N,y0,abstol,reltol,'PI');
[tnList4,ynList4,hList4,rList4,nfun4] = ImplicitTrapezoidAdaptiveStep(...
    @VanDerPol,@jacVanDerPol,tspan,N,y0,abstol,reltol,'PI');

% Plotting
figure(1)
lw = 2;
subplot(2,1,1)
hold on
    plot(tnList1(1,:),ynList1(1,:),'-.','LineWidth',lw)
    plot(tnList4(1,:),ynList4(1,:),':','LineWidth',lw)
    leg = legend('Explicit Euler','Trapezoidal');
    set(leg,'FontSize',14);
    ylabel('x','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

subplot(2,1,2)
hold on
    plot(tnList1(1,:),ynList1(2,:),'-.','LineWidth',lw)
    plot(tnList4(1,:),ynList4(2,:),':','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4');
    set(leg,'FontSize',14);
    ylabel('y','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

figure(2)
lw = 2;
subplot(2,1,1)
hold on
    plot(tnList1(2:end),hList1,'-.','LineWidth',lw)
    plot(tnList2(2:end),hList2,'-.','LineWidth',lw)
    plot(tnList3(2:end),hList3,':','LineWidth',lw)
    plot(tnList4(2:end),hList4,'--','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4');
    set(leg,'FontSize',14);
    ylabel('h','FontSize',14,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',14,'FontWeight','bold')
    grid on
hold off

subplot(2,1,2)
ve = ones(1,length(tnList1)-1);
hold on
    plot(tnList1(2:end),rList1,'-.','LineWidth',lw)
    plot(tnList2(2:end),rList2,'-.','LineWidth',lw)
    plot(tnList3(2:end),rList3,':','LineWidth',lw)
    plot(tnList4(2:end),rList4,'--','LineWidth',lw)
    plot(tnList1(2:end),ve,'LineWidth',3)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4','MaxTolerance');
    set(leg,'FontSize',14);
    ylim([0 1.3])
    ylabel('r','FontSize',14,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',14,'FontWeight','bold')
    grid on
hold off
