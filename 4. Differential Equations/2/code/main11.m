%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exercise 1.1: Newton?s method for solving nonlinear BVPs  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subscripts
addpath('functions');
%% Ex 1.1.3 BVP 
epsilon = 0.05;
twoode = @(t,x,epsilon) [x(2); -x(1)*(x(2)-1)/epsilon];
twobc = @(xa,xb,epsilon) [xa(1)+1; xb(1)-1.5];
options=bvpset('reltol',1e-2,'abstol',1e-2);
solinit = bvpinit(linspace(0,1,10),[-1 1.5]);
sol = bvp4c(twoode,twobc,solinit,options,epsilon);

figure(1)
hold on 
    plot(sol.yp(2,:))
    grid on
hold off

%% Ex 1.1.4: BVP Finite Differences Method for non-linear equations

close all
%System parameters
m = 400;                            %number of points on grid
h = 1/(m);
x = linspace(0,1,m+1);                %x gird
epsilon = 0.001;     
alpha = -1;                         %u(0)=alpha
beta = 1.5;                         %u(1)=beta
params = [alpha beta epsilon];

%Analytic guess of solution
a = 0;
b = 1;
w0 = 1/2*(a-b+beta-alpha);
x_hat = 1/2*(a+b-alpha-beta);
u_guess = x-x_hat+w0*tanh(w0*(x-x_hat)/(2*epsilon));

%Newtons method
U = newtonODE(@functions11,u_guess',h,10e-8,epsilon);

% Plotting
lw = 2;  
plot(x,U,'LineWidth',lw)
hold on    
    plot(x,u_guess,'--','LineWidth',lw+1)
%     plot(x(2:n-1),x(2:n-1)+alpha-a,'-.','LineWidth',lw)
%     plot(x(2:n-1),x(2:n-1)+beta-b,'-.','LineWidth',lw)
    leg = legend('Estimated solution','Analytical guess');
    set(leg,'FontSize',14);
    ep = num2str(epsilon);
    title(['BVP for epsilon = ' ep],'FontSize',14);
    ylabel('u(x)','FontSize',14)
    xlabel('x','FontSize',14)
    grid on
hold off

