%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Exercise 1.2: Shooting Methods  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add paths
addpath('functions');

%% Ex 1.2.1 Shooting methods through Newtons root finding

close all
%System paramters
alpha = -1;          % u(0)
beta = 1.5;         % u(1)
sigma0 = 1.001;       % u'(0) guess
epsilon = 0.001;     
sigma1 = sigma0;
%Solver parameters
options = odeset('RelTol',1e-6);
tspan = [0 1];
%Newton method paramters
tol = 10e-5;
err = 1000;
itmax = 1000;
it = 0;
%Newtons method begins
while(err>tol && it<itmax)
    sigma = sigma1;
    Y0 = [alpha sigma 0 1];
    [t,y] = ode45(@(x,u) func12Newton(x,u,epsilon), tspan, Y0, options);
    G = y(end,1)-beta;
    dG = y(end,3);
    sigma1 = sigma - G/dG;
    err = norm(y(end,1)-beta);
    it = it+1;
end

figure(1)
plot(t,y(:,1),'LineWidth',lw)
hold on
    ep = num2str(epsilon);
    title(['Solution for epsilon = ' ep],'FontSize',14);
    ylabel('U(x)','FontSize',14)
    xlabel('x','FontSize',14)
    grid
hold off


%% Ex 1.2.2 Sensivity Analysis 

close all
%Memory allocation
n = 100;
uList = zeros(n,1);                 % u(1) for different sigmas
wList = zeros(n,1);                 % du(1)/dsigma for different sigmas
%System parameters
alpa = -1;                          % u(0)
beta = 1.5;                         % u(1)
sigmaList = linspace(1.1,10,n);     % u'(0) guess list
epsilon = 0.001;
%Solver parametes
options = odeset('RelTol',1e-6);
tspan = [0 1];

for i = 1:n
    Y0 = [alpa sigmaList(i) 0 1];
    [t,y] = ode45(@(x,u) func12Newton(x,u,epsilon), tspan, Y0, options);
    uList(i) = y(end,1);
    wList(i) = y(end,3);
end

%Plotting
lw = 2;
figure(1)
plot(sigmaList,wList,'LineWidth',lw)
hold on
    ep = num2str(epsilon);
    title(['Sensivity Analysis for epsilon = ' ep],'FontSize',14);
    ylabel('S_\sigma','FontSize',14)
    xlabel('u´(0)=\sigma','FontSize',14)
    grid
hold off

figure(2)
plot(sigmaList,uList,'LineWidth',lw)
hold on
    grid
    title(['u(1,\sigma) vs Sigma for epsilon = ' ep],'FontSize',14);
    ylabel('u(1,\sigma)','FontSize',14)
    xlabel('u´(0)=\sigma','FontSize',14)
hold off

%% Ex 1.2.1 With psudo-secant method

close all
alpa = -1;      % u(0)
beta = 1.5;     % u(1)
sigma = 1.1;    % u'(0) guess
epsilon = 0.001;
Y0 = [alpha sigma];
options = odeset('RelTol',1e-6);
tspan = [0 1];
tol = 10e-6;

error = 1000;
c = [1.1 1.3];
while(error > tol)    
    [t,y] = ode45(@(x,u) func12(x,u,epsilon), tspan, [alpha c(1)], options);
    c1 = y(end,1)-beta;
    [t,y] = ode45(@(x,u) func12(x,u,epsilon), tspan, [alpha c(2)], options);
    c2 = y(end,1)-beta;
    l = 0.5*(c(1)+c(2));
    [t,y] = ode45(@(x,u) func12(x,u,epsilon), tspan, [alpha l], options);
    c3 = y(end,1)-beta;
    if(c3 < 0)
        c(1) = l;
    else
        c(2) = l;
    end
    error = norm(c1-c2);
end
[t,y] = ode45(@(x,u) func12(x,u,epsilon), tspan, [alpha c(2)], options);
plot(t,y(:,1))


