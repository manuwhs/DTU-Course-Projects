%% Question 1: mu = 3
close all; clear all;

% Paramters
tspan = [0 100];
N = 100000;
y0 = [1;1];
method = 'DOPRI54';
solver = ERKSolverErrorEstimationParameters(method);
h = (tspan(2)-tspan(1))/N

% Solvers
[tnList1,ynList1] = ExplicitEulers(@VanDerPol,tspan,N,y0);
[tnList2,ynList2] = ImplicitEulers(@VanDerPol,@jacVanDerPol,tspan,N,y0);
[tnList3,ynList3] = ImplicitTrapezoid(@VanDerPol,@jacVanDerPol,tspan,N,y0);
[tnList4,ynList4] = RK4(@VanDerPol,tspan,N,y0);
[tnList5,ynList5,Eout5] = ...
    ERKSolverErrorEstimation(@VanDerPol,tspan,y0,h,solver);
ynList5 = ynList5';

% Double of steps for making error estimation
[tnList11,ynList11] = ExplicitEulers(@VanDerPol,tspan,N*2,y0);
[tnList22,ynList22] = ImplicitEulers(@VanDerPol,@jacVanDerPol,tspan,N*2,y0);
[tnList33,ynList33] = ImplicitTrapezoid(@VanDerPol,@jacVanDerPol,tspan,N*2,y0);
[tnList44,ynList44] = RK4(@VanDerPol,tspan,N*2,y0);
[tnList55,ynList55,Eout55] = ...
    ERKSolverErrorEstimation(@VanDerPol,tspan,y0,h/2,solver);
ynList55 = ynList55';

len = length(ynList1);
LTEest1 = zeros(1,len);
LTEest2 = zeros(1,len);
LTEest3 = zeros(1,len);
LTEest4 = zeros(1,len);
LTEest5 = zeros(1,len);
LTEembed5 = zeros(1,len);
for i = 2:len
    LTEest1(i) = 2*norm(ynList1(:,i)-ynList11(:,i+i-1));
    LTEest2(i) = 2*norm(ynList2(:,i)-ynList22(:,i+i-1));
    LTEest3(i) = (4/3)*norm(ynList3(:,i)-ynList33(:,i+i-1));
    LTEest4(i) = (16/15)*norm(ynList4(:,i)-ynList44(:,i+i-1));
    LTEest5(i) = (32/31)*norm(ynList5(:,i)-ynList55(:,i+i-1));
    LTEembed5(i) = norm(Eout5(2,:));
end
%% Plotting question 1
lw = 2;
figure(1)
hold on
    title('Van Der Pol mu = 100');
    plot(ynList1(1,:),ynList1(2,:),'-.','LineWidth',lw)
    plot(ynList2(1,:),ynList2(2,:),':','LineWidth',lw)
    plot(ynList3(1,:),ynList3(2,:),'--','LineWidth',lw)
    plot(ynList4(1,:),ynList4(2,:),'--','LineWidth',lw+1)
    plot(ynList5(1,:),ynList5(2,:),':','LineWidth',lw+2)
    leg = legend('Explicit Euler','Implicit euler',...
        'Trapezoidal','RK4','DOPRI54');
    set(leg,'FontSize',12);
    ylabel('y(t)','FontSize',12,'FontWeight','bold')
    xlabel('x(t)','FontSize',12,'FontWeight','bold')
    grid on
hold off

figure(2)
semilogy(tnList1,LTEest1,'LineWidth',lw)
hold on
    title('Estimated error for mu = 100');
    semilogy(tnList1,LTEest2,'LineWidth',lw)
    semilogy(tnList1,LTEest3,'LineWidth',lw)
    semilogy(tnList1,LTEest4,'LineWidth',lw)
    semilogy(tnList1,LTEest5,'LineWidth',lw)
    leg = legend('Explicit Euler','Implicit euler',...
        'Trapezoidal','RK4','DOPRI54');
    set(leg,'FontSize',12);
    ylabel('log(LTEest)','FontSize',12,'FontWeight','bold')
    xlabel('time(seconds)','FontSize',12,'FontWeight','bold')
    grid on
hold off
%%
NList = 10:1000;
len = length(NList);
hList = zeros(1,len);
LTEest1 = zeros(1,len);
LTEest2 = zeros(1,len);
LTEest3 = zeros(1,len);
LTEest4 = zeros(1,len);
LTEest5 = zeros(1,len);
LTEembed5 = zeros(1,len);
for i = 1:len
    h = (tspan(2)-tspan(1))/NList(i);
    hList(i) = h;
    [tnList1,ynList1] = ExplicitEulers(@VanDerPol,tspan,NList(i),y0);
    [tnList11,ynList11] = ExplicitEulers(@VanDerPol,tspan,NList(i),y0);
    [tnList2,ynList2] = ImplicitEulers(@VanDerPol,@jacVanDerPol,tspan,NList(i),y0);
    [tnList22,ynList22] = ImplicitEulers(@VanDerPol,@jacVanDerPol,tspan,NList(i),y0);
    [tnList3,ynList3] = ExplicitTrapezoid(@VanDerPol,tspan,NList(i),y0);
    [tnList33,ynList33] = ExplicitTrapezoid(@VanDerPol,tspan,NList(i),y0);
    [tnList4,ynList4] = RK4(@VanDerPol,tspan,NList(i),y0);
    [tnList44,ynList44] = RK4(@VanDerPol,tspan,NNList(i),y0);
    [tnList5,ynList5,Eout1] = ...
    ERKSolverErrorEstimation(@VanDerPol,tspan,y0,h,solver,lambda);
    [tnList55,ynList55,Eout2] = ...
    ERKSolverErrorEstimation(@VanDerPol,tspan,y0,h,solver,lambda);

    LTEest1(i) = 2*norm(ynList1(:,2)-ynList11(:,3));
    LTEest2(i) = 2*norm(ynList2(:,2)-ynList22(:,3));
    LTEest3(i) = (4/3)*norm(ynList3(:,2)-ynList33(:,3));
    LTEest4(i) = (16/15)*norm(ynList4(:,2)-ynList44(:,3));
    LTEest5(i) = (32/31)*norm(ynList5(:,2)-ynList55(:,3));
    LTEembed5(i) = norm(Eout(2,:));
end