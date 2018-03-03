%% Question 3: IVP1 global error until time t=10
close all; clear all;

% Paramters
tspan = [0 10];
N = 40;
y0 = 1;
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
lambda = -1;
h = (tspan(2)-tspan(1))/N;

% Solvers
[tnList,ynList] = ExplicitEulers(@func,tspan,N,y0);
[tnList1,ynList1] = ImplicitEulers(@func,@Jacob,tspan,N,y0);
[tnList2,ynList2] = ImplicitTrapezoid(@func,@Jacob,tspan,N,y0);
[tnList3,ynList3] = RK4(@func,tspan,N,y0);
[tnList4,ynList4,Eout] = ...
    ERKSolverErrorEstimation(@func,tspan,y0,h,solver,lambda);

% Plotting all solvers in log scale
y = exp(lambda*tnList)*y0;
y1 = exp(lambda*tnList4)*y0;
lw = 2;
semilogy(tnList,abs(y-ynList),'LineWidth',lw)
hold on
    title('\textbf{Global error over time for} $\dot{x}(t) = -x(t)$',...
        'interpreter','latex');
    semilogy(tnList1,abs(y-ynList1),'LineWidth',lw)
    semilogy(tnList2,abs(y-ynList2),'LineWidth',lw)
    semilogy(tnList3,abs(y-ynList3),'LineWidth',lw)
    semilogy(tnList4,abs(y1-ynList4),'LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4','DOPRI54');
    set(leg,'FontSize',14);
    ylabel('log(error)','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid
hold off

% Plotting Explicit vs Implicit Euler with and without log scale
figure(2)
subplot(2,1,1);  title('Global error: Explicit and Implicit Euler');
semilogy(tnList,abs(y-ynList),'LineWidth',lw)
hold on
    semilogy(tnList1,abs(y-ynList1),'LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('log(Global error)','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off
subplot(2,1,2); title('Global error: without loagarithmic scale');
hold on
    plot(tnList,abs(y-ynList),'LineWidth',1)
    plot(tnList1,abs(y-ynList1),'LineWidth',1)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('Global error','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

figure(3)
hold on
    title('\textbf{Solution for} $\dot{x}(t) = -x(t)$',...
        'interpreter','latex');
    plot(tnList,ynList,'--','LineWidth',lw)
    plot(tnList,ynList1,':','LineWidth',lw)
    plot(tnList,ynList2,'-.','LineWidth',lw)
    plot(tnList,ynList3,'--','LineWidth',lw)
    plot(tnList,ynList4,':','LineWidth',4)
    plot(tnList,y,':','LineWidth',5)
    leg = legend('Explicit Euler','Implicit Euler','Trapezoid',...
        'RK4','DOPRI54','Real');
    set(leg,'FontSize',14);
    ylabel('x(t)','FontSize',12,'FontWeight','bold')
    xlabel('t(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

%% Question 3: IVP2 global error until time t=10
close all; clear all;
tspan = [0 10];
N = 40;
Y0 = [0 ;1];
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
h = (tspan(2)-tspan(1))/N;

% Solvers
[tnList,ynList] = ExplicitEulers(@IVP2,tspan,N,Y0);
[tnList1,ynList1] = ImplicitEulers(@IVP2,@jacIVP2,tspan,N,Y0);
[tnList2,ynList2] = ImplicitTrapezoid(@IVP2,@jacIVP2,tspan,N,Y0);
[tnList3,ynList3] = RK4(@IVP2,tspan,N,Y0);
[tnList4,ynList4,Eout] = ...
    ERKSolverErrorEstimation(@IVP2,tspan,Y0,h,solver);
ynList4 = ynList4';

% Analytical solution
A = [0 1;-1 0];
[m,n] = size(Y0);
timeLen = length(tnList);
Y = zeros(m,timeLen);
for i = 1:timeLen
    Y(:,i) = expm(A.*tnList(i))*Y0;
end

% Plotting
figure(1)
lw = 2;
semilogy(tnList,globError(Y,ynList),'LineWidth',lw)
hold on
    title('\textbf{Global error over time for} $\ddot{x}(t) = -x(t)$',...
        'interpreter','latex');
    semilogy(tnList,globError(Y,ynList1),'LineWidth',lw)
    semilogy(tnList,globError(Y,ynList2),'LineWidth',lw)
    semilogy(tnList,globError(Y,ynList3),'LineWidth',lw)
    semilogy(tnList,globError(Y,ynList4),'LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4','DOPRI54');
    set(leg,'FontSize',14);
    ylabel('log(Global_error)','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid
hold off

% Plotting Explicit vs Implicit Euler in with and without log scale
figure(2)
subplot(2,1,1);  title('Global error: Explicit and Implicit Euler');
semilogy(tnList,globError(Y,ynList),'LineWidth',1)
hold on
    semilogy(tnList,globError(Y,ynList1),'LineWidth',1)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('log(Global error)','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off
subplot(2,1,2); title('Global error: without logarithmic scale');
hold on
    plot(tnList,globError(Y,ynList),'LineWidth',1)
    plot(tnList,globError(Y,ynList1),'LineWidth',1)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('Global error','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

figure(3)
hold on
    title('\textbf{Solution for} $\ddot{x}(t) = -x(t)$',...
        'interpreter','latex');
    plot(tnList,ynList(1,:),'--','LineWidth',lw)
    plot(tnList,ynList1(1,:),':','LineWidth',lw)
    plot(tnList,ynList2(1,:),'-.','LineWidth',lw)
    plot(tnList,ynList3(1,:),'--','LineWidth',lw)
    plot(tnList,ynList4(1,:),':','LineWidth',4)
    plot(tnList,Y(1,:),':','LineWidth',5)
    leg = legend('Explicit Euler','Implicit Euler','Trapezoid',...
        'RK4','DOPRI54','Real');
    set(leg,'FontSize',14);
    ylabel('x(t)','FontSize',12,'FontWeight','bold')
    xlabel('t(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off

%% Question 4: LTE error for differnet values of h IVP2
close all; clear all;
tspan = [0 10];
Y0 = [0 ;1];
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
A = [0 1;-1 0];

NList = 10:1000;
len = length(NList);
hList = zeros(1,len);
LTE = zeros(1,len);
LTE1 = zeros(1,len);
LTE2 = zeros(1,len);
LTE3 = zeros(1,len);
LTE4 = zeros(1,len);
for i = 1:len
    h = (tspan(2)-tspan(1))/NList(i);
    hList(i) = h;
    [tnList,ynList] = ExplicitEulers(@IVP2,tspan,NList(i),Y0);
    [tnList1,ynList1] = ImplicitEulers(@IVP2,@jacIVP2,tspan,NList(i),Y0);
    [tnList2,ynList2] = ImplicitTrapezoid(@IVP2,@jacIVP2,tspan,NList(i),Y0);
    [tnList3,ynList3] = RK4(@IVP2,tspan,NList(i),Y0);
    [tnList4,ynList4,Eout] = ...
        ERKSolverErrorEstimation(@IVP2,tspan,Y0,h,solver);
    ynList4 = ynList4';
    Y = expm(A.*tnList(2))*Y0;

    LTE(i) = norm(ynList(:,2)-Y);
    LTE1(i) = norm(ynList1(:,2)-Y);
    LTE2(i) = norm(ynList2(:,2)-Y);
    LTE3(i) = norm(ynList3(:,2)-Y);
    LTE4(i) = norm(ynList4(:,2)-Y);
end
%%
m1= (log(LTE(100))-log(LTE(200)))/(log(hList(100))-log(hList(200)));
m2= (log(LTE1(100))-log(LTE1(200)))/(log(hList(100))-log(hList(200)));
m3= (log(LTE2(100))-log(LTE2(200)))/(log(hList(100))-log(hList(200)));
m4= (log(LTE3(100))-log(LTE3(200)))/(log(hList(100))-log(hList(200)));
m5= (log(LTE4(100))-log(LTE4(200)))/(log(hList(100))-log(hList(200)));
%%
close all;
lw = 2;
figure(1)
loglog(hList,LTE,'LineWidth',lw)
hold on
    title('\textbf{LTE vs step size for} $\ddot{x}(t) = -x(t)$',...
        'interpreter','latex','FontSize',14);
    
    loglog(hList,LTE1,'--','LineWidth',lw)
    loglog(hList,LTE2,'LineWidth',lw)
    loglog(hList,LTE3,'LineWidth',lw)
    loglog(hList,LTE4,'LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4','DOPRI54');
    set(leg,'FontSize',14);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off

% Plotting Explicit vs Implicit Euler in with and without log scale
figure(2)
subplot(2,1,1);  title('LTE: Explicit and Implicit Euler');
loglog(hList,LTE,'LineWidth',lw)
hold on
    loglog(hList,LTE1,'--','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off
subplot(2,1,2); title('LTE: without logarithmic scale');
hold on
    plot(hList,LTE,'LineWidth',lw)
    plot(hList,LTE1,'--','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('LTE error','FontSize',12,'FontWeight','bold')
    xlabel('h','FontSize',12,'FontWeight','bold')
    grid on
hold off

%% Question 4: LTE error for differnet values of h IVP1
close all; clear all;
tspan = [0 10];
y0 = 1;
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
lambda = -1;

NList = 10:1000;
len = length(NList);
hList = zeros(1,len);
LTE = zeros(1,len);
LTE1 = zeros(1,len);
LTE2 = zeros(1,len);
LTE3 = zeros(1,len);
LTE4 = zeros(1,len);
for i = 1:len
    h = (tspan(2)-tspan(1))/NList(i);
    hList(i) = h;
    [tnList,ynList] = ExplicitEulers(@func,tspan,NList(i),y0);
    [tnList1,ynList1] = ImplicitEulers(@func,@Jacob,tspan,NList(i),y0);
    [tnList2,ynList2] = ImplicitTrapezoid(@func,@Jacob,tspan,NList(i),y0);
    [tnList3,ynList3] = RK4(@func,tspan,NList(i),y0);
    [tnList4,ynList4,Eout] = ...
    ERKSolverErrorEstimation(@func,tspan,y0,h,solver,lambda);
    ynList4 = ynList4';

    y = exp(lambda*tnList(2))*y0;

    LTE(i) = norm(ynList(:,2)-y);
    LTE1(i) = norm(ynList1(:,2)-y);
    LTE2(i) = norm(ynList2(:,2)-y);
    LTE3(i) = norm(ynList3(:,2)-y);
    LTE4(i) = norm(ynList4(:,2)-y);
end
%%
m1= (log(LTE(100))-log(LTE(200)))/(log(hList(100))-log(hList(200)));
m2= (log(LTE1(100))-log(LTE1(200)))/(log(hList(100))-log(hList(200)));
m3= (log(LTE2(100))-log(LTE2(200)))/(log(hList(100))-log(hList(200)));
m4= (log(LTE3(100))-log(LTE3(200)))/(log(hList(100))-log(hList(200)));
m5= (log(LTE4(100))-log(LTE4(200)))/(log(hList(100))-log(hList(200)));
%%
close all;
lw = 2;
loglog(hList,LTE,'LineWidth',lw)
hold on
    title('\textbf{LTE vs step size for} $\dot{x}(t) = -x(t)$',...
        'interpreter','latex','FontSize',14);
    
    loglog(hList,LTE1,'--','LineWidth',lw)
    loglog(hList,LTE2,'LineWidth',lw)
    loglog(hList,LTE3,'LineWidth',lw)
    loglog(hList,LTE4,'LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler',...
        'Trapezoidal','RK4','DOPRI54');
    set(leg,'FontSize',14);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off

% Plotting Explicit vs Implicit Euler in with and without log scale
figure()
subplot(2,1,1);  
loglog(hList,LTE,'LineWidth',lw)
hold on
    title('LTE: Explicit and Implicit Euler')
    loglog(hList,LTE1,'--','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off
subplot(2,1,2); title('LTE: without logarithmic scale');
hold on
    plot(hList,LTE,'LineWidth',lw)
    plot(hList,LTE1,'--','LineWidth',lw)
    leg = legend('Explicit Euler','Implicit Euler');
    set(leg,'FontSize',14);
    ylabel('LTE error','FontSize',12,'FontWeight','bold')
    xlabel('h','FontSize',12,'FontWeight','bold')
    grid on
hold off

% Calculate slopes for order checking
m1= (log(LTE(10))-log(LTE(30)))/(log(hList(10))-log(hList(30)));
m2= (log(LTE1(10))-log(LTE1(30)))/(log(hList(10))-log(hList(30)));
m3= (log(LTE2(10))-log(LTE2(30)))/(log(hList(10))-log(hList(30)));
m4= (log(LTE3(10))-log(LTE3(30)))/(log(hList(10))-log(hList(30)));
m5= (log(LTE4(10))-log(LTE4(30)))/(log(hList(10))-log(hList(30)));
 %% Question 5: error estimation for IVP1
close all; clear all;
tspan = [0 10];
y0 = 1;
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
lambda = -1;

NList = 10:1000;
len = length(NList);
hList = zeros(1,len);
LTEest1 = zeros(1,len);
LTEreal1 = zeros(1,len);
LTEest2 = zeros(1,len);
LTEreal2 = zeros(1,len);
LTEest3 = zeros(1,len);
LTEreal3 = zeros(1,len);
LTEest4 = zeros(1,len);
LTEreal4 = zeros(1,len);
LTEest5 = zeros(1,len);
LTEreal5 = zeros(1,len);
for i = 1:len
    h = (tspan(2)-tspan(1))/NList(i);
    hList(i) = h;
    [tnList1,ynList1] = ExplicitEulers(@func,tspan,NList(i),y0);
    [tnList11,ynList11] = ExplicitEulers(@func,tspan,2*NList(i),y0);
    [tnList2,ynList2] = ImplicitEulers(@func,@Jacob,tspan,NList(i),y0);
    [tnList22,ynList22] = ImplicitEulers(@func,@Jacob,tspan,2*NList(i),y0);
    [tnList3,ynList3] = ImplicitTrapezoid(@func,@Jacob,tspan,NList(i),y0);
    [tnList33,ynList33] = ImplicitTrapezoid(@func,@Jacob,tspan,2*NList(i),y0);
    [tnList4,ynList4] = RK4(@func,tspan,NList(i),y0);
    [tnList44,ynList44] = RK4(@func,tspan,2*NList(i),y0);
  
    y = exp(lambda*tnList1(2))*y0;
    LTEest1(i) = 2*norm(ynList1(2)-ynList11(3));
    LTEreal1(i) = norm(ynList1(2)-y);
    LTEest2(i) = 2*norm(ynList2(2)-ynList22(3));
    LTEreal2(i) = norm(ynList2(2)-y);
    LTEest3(i) = (4/3)*norm(ynList3(2)-ynList33(3));
    LTEreal3(i) = norm(ynList3(2)-y);
    LTEest4(i) = (16/15)*norm(ynList4(2)-ynList44(3));
    LTEreal4(i) = norm(ynList4(2)-y);
end
%%
figure(1)
loglog(hList,LTEest1,':','LineWidth',3)
hold on
    loglog(hList,LTEreal1,'--','LineWidth',2)
    loglog(hList,LTEest2,':','LineWidth',2)
    loglog(hList,LTEreal2,'--','LineWidth',2)
    loglog(hList,LTEest3,':','LineWidth',2)
    loglog(hList,LTEreal3,'--','LineWidth',2)
    loglog(hList,LTEest4,':','LineWidth',2)
    loglog(hList,LTEreal4,'--','LineWidth',3)
    leg = legend('Error estimation ExEuler','Real error ExEuler',...
        'Error estimation Trapezoid','Real error Trapezoid',...
        'Error estimation RK4','Real error RK4');
    set(leg,'FontSize',12);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off
    
%% Question 5: error estimation for IVP2
close all; clear all;
tspan = [0 10];
Y0 = [0 ;1];
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
A = [0 1;-1 0];

NList = 10:1000;
len = length(NList);
hList = zeros(1,len);
LTEest1 = zeros(1,len);
LTEreal1 = zeros(1,len);
LTEest2 = zeros(1,len);
LTEreal2 = zeros(1,len);
LTEest3 = zeros(1,len);
LTEreal3 = zeros(1,len);
LTEest4 = zeros(1,len);
LTEreal4 = zeros(1,len);
LTEest5 = zeros(1,len);
LTEreal5 = zeros(1,len);
LTEembed5 = zeros(1,len);
for i = 1:len
    h = (tspan(2)-tspan(1))/NList(i);
    hList(i) = h;
    [tnList1,ynList1] = ExplicitEulers(@IVP2,tspan,NList(i),Y0);
    [tnList11,ynList11] = ExplicitEulers(@IVP2,tspan,2*NList(i),Y0);
    [tnList2,ynList2] = ImplicitEulers(@IVP2,@jacIVP2,tspan,NList(i),Y0);
    [tnList22,ynList22] = ImplicitEulers(@IVP2,@jacIVP2,tspan,2*NList(i),Y0);
    [tnList3,ynList3] = ImplicitTrapezoid(@IVP2,@jacIVP2,tspan,NList(i),Y0);
    [tnList33,ynList33] = ImplicitTrapezoid(@IVP2,@jacIVP2,tspan,2*NList(i),Y0);
    [tnList4,ynList4] = RK4(@IVP2,tspan,NList(i),Y0);
    [tnList44,ynList44] = RK4(@IVP2,tspan,2*NList(i),Y0);
    Y = expm(A.*tnList1(2))*Y0;
    [tnList5,ynList5,Eout] = ...
        ERKSolverErrorEstimation(@IVP2,tspan,Y0,h,solver);
    ynList5 = ynList5';
    [tnList55,ynList55,Eout1] = ...
        ERKSolverErrorEstimation(@IVP2,tspan,Y0,h/2,solver);
    ynList55 = ynList55';

    LTEest1(i) = 2*norm(ynList1(:,2)-ynList11(:,3));
    LTEreal1(i) = norm(ynList1(:,2)-Y);
    LTEest2(i) = 2*norm(ynList2(:,2)-ynList22(:,3));
    LTEreal2(i) = norm(ynList2(:,2)-Y);
    LTEest3(i) = (4/3)*norm(ynList3(:,2)-ynList33(:,3));
    LTEreal3(i) = norm(ynList3(:,2)-Y);
    LTEest4(i) = (16/15)*norm(ynList4(:,2)-ynList44(:,3));
    LTEreal4(i) = norm(ynList4(:,2)-Y);
    LTEest5(i) = (32/31)*norm(ynList5(:,2)-ynList55(:,3));
    LTEreal5(i) = norm(ynList5(:,2)-Y);
    LTEembed5(i) = norm(Eout(2,:));
end
%%
close all;
figure(1)
loglog(hList,LTEest1,':','LineWidth',3)
hold on
    title('Error estimaton by step-doubling')
    loglog(hList,LTEreal1,'--','LineWidth',2)
    loglog(hList,LTEest2,':','LineWidth',2)
    loglog(hList,LTEreal2,'--','LineWidth',2)
    loglog(hList,LTEest3,':','LineWidth',2)
    loglog(hList,LTEreal3,'--','LineWidth',2)
    loglog(hList,LTEest4,':','LineWidth',2)
    loglog(hList,LTEreal4,'--','LineWidth',3)
    leg = legend('Error estimation ExEuler','Real error ExEuler',...
        'Error estimation ImpEuler','Real error ImEuler',...
        'Error estimation Trapezoid','Real error Trapezoid',...
        'Error estimation RK4','Real error RK4');
    set(leg,'FontSize',12);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off

figure(2)
loglog(hList,LTEreal5,':','LineWidth',4)
hold on
    title('DOPRI error estimation')
    loglog(hList,LTEembed5,'--','LineWidth',2)
    loglog(hList,LTEest5,'-.','LineWidth',2)
    leg = legend('Real error','Embbeded estimation',...
        'Step-doubling Estimation');
    set(leg,'FontSize',12);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off

figure(3)
loglog(hList,LTEest1,':','LineWidth',2)
hold on
    title('Explicit and Implicit Euler error estimation')
    loglog(hList,LTEreal1,'--','LineWidth',2)
    loglog(hList,LTEest2,':','LineWidth',2)
    loglog(hList,LTEreal2,'--','LineWidth',2)
    leg = legend('Error estimation ExEuler','Real error ExEuler',...
        'Error estimation ImpEuler','Real error ImEuler');
    set(leg,'FontSize',12);
    ylabel('log(LTE)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off