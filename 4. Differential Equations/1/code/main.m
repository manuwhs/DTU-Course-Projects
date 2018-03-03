%% Explicit eulers
close all; clear all;

tspan = 10;
Nsteps = 20;
lw = 3;
[tnList_0,ynList_0] = ExplicitEulers(@func,tspan,Nsteps*2,1);
[tnList_1,ynList_1] = ExplicitTrapezoid(@func,tspan,Nsteps,1);
[tnList_2,ynList_2] = Midpoint(@func,tspan,Nsteps,1);
[tnList_4,ynList_4] = RK4(@func,tspan,Nsteps,1);

%% Error
figure()
subplot(2,1,1);  title('Error Explicit Euler');
Error = ynList_0 - exp(1*tnList_0)*1;
plot(tnList_0, Error, 'color',rand(1,3), 'LineWidth',lw)
subplot(2,1,2); title('Error RK4');
Error = tnList_4 - exp(1*tnList_4)*1;
plot(tnList_4, Error, 'color',rand(1,3), 'LineWidth',lw)
figure()



hold on
    plot(tnList_0, ynList_0, 'color',rand(1,3), 'LineWidth',lw)
    plot(tnList_1, ynList_1, 'color',rand(1,3),'LineWidth',lw)
    plot(tnList_2, ynList_2, 'color',rand(1,3),'LineWidth',lw)
    
    plot(tnList_4, ynList_4, 'color',rand(1,3),'LineWidth',lw)
    plot(tnList_0, exp(1*tnList_0)*1, 'color',rand(1,3),'LineWidth',lw)  %e^(lambda*x)*x0
    legend('Euler', 'Trapezoid','Midpint','RK4','real')
hold off

N = 10;
alphas = (1:N-1)/N;

figure()
hold on
for i = 1:N-1
    [tnList,ynList] = RK2(@func,tspan,Nsteps,1, 0.4);
    plot(tnList, ynList, 'color',rand(1,3), 'LineWidth',lw, 'DisplayName',num2str(alphas(i)))
end 
hold off

%% Checking 2D work
tspan = 50;
Y0 = [2 ; 2];
[tnList_X,ynList_X] = RK4(@DepPrey,tspan,Nsteps*1000,Y0);  % DepPrey

figure()

subplot(2,1,1);  
plot(ynList_X(1,:), ynList_X(2,:), 'color',rand(1,3), 'LineWidth',lw)

mu = 3;
tspan = [0 5*mu];
Y0 = [2 ; 0];
[tnList_X,ynList_X] = RK4(@VanDelPol,tspan,Nsteps*1000,Y0);  

subplot(2,1,2);  
plot(ynList_X(1,:), ynList_X(2,:), 'color',rand(1,3), 'LineWidth',lw)

%% Checking tstar and tend work
tspan = 10;
Nsteps = 50;
[tnList_0,ynList_0] = ExplicitTrapezoid(@func,tspan,Nsteps,1);
[tnList_1,ynList_1] = ExplicitTrapezoid(@func,[1  tspan+1],Nsteps,1);
[tnList_2,ynList_2] = ExplicitTrapezoid(@func,[2  tspan+2],Nsteps,1);

figure()
hold on
    plot(tnList_0, ynList_0, 'color',rand(1,3), 'LineWidth',lw)
    plot(tnList_1, ynList_1, 'color',rand(1,3),'LineWidth',lw)
    plot(tnList_2, ynList_2, 'color',rand(1,3),'LineWidth',lw)
hold off

%% Implicit Euler 1
tspan = 5;
Nsteps = 100;
[tnList_0,ynList_0] = ExplicitEulers(@func,tspan,Nsteps,1);
[tnList_1,ynList_1] = ImplicitEulers1(@func,tspan,Nsteps,1, 100);


figure()
hold on
    plot(tnList_0, ynList_0, 'color',rand(1,3), 'LineWidth',lw)
    plot(tnList_1, ynList_1, 'color',rand(1,3), 'LineWidth',lw)
    plot(tnList_0, exp(1*tnList_0)*1, 'color',rand(1,3),'LineWidth',lw)  %e^(lambda*x)*x0

    legend('Explicit','Implicit','real')
hold off
% line width 3
% font size 16

%% Impicit Eulers