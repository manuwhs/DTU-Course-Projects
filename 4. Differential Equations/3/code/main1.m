%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Exercise 1: Parabolic equations   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subscripts
addpath('functions');
%% Forward time Central Space
%model paramters
kappa = 0.1;

% Space grid
m = 100;
h = 2/(m-1);
x = linspace(-1,1,m);
%% Forward time Central Space
%model paramters
theta = 0;
% to fulfil the stability criterion k > 0.5*h^2/kappa
k = 0.4*h^2/kappa;      % good k
% k = h^2/(6*kappa);      % extremely good k   
% k = 0.6*h^2/kappa;	% wrong k

%Solving the system
steps = 4000;
Un = solver(m,k,h,theta,x,steps);
Uf = realValue(m,x,steps*k);
U0 = realValue(m,x,0);

%Convergence
k_ex = @(h,kappa) h^2/(6*kappa);
k_good = @(h,kappa) 0.4*h^2/kappa;
k_wrong = @(h,kappa) 0.6*h^2/kappa;
m_List = 6:5:200;
[LTE_List,h_List,k_List] = LTEconvergence(m_List,kappa,...
    theta,k_ex);
[LTE_List1,h_List1,k_List1] = LTEconvergence(m_List,kappa,...
    theta,k_good);

%%%%% Plotting

%Solution after time
figure(1)
lw = 2;
plot(U0,'-.','LineWidth',lw)
hold on
    plot(Un,'LineWidth',lw)
    plot(Uf,'-.','LineWidth',lw)
    leg = legend('U_0','U_f estimated','U_f real');
    set(leg,'FontSize',14);
    ste = num2str(steps);
    kstring = num2str(k);
    title(['Solution for n = ' ste ' and k = ' kstring],'FontSize',14);
    ylabel('u(x)','FontSize',14)
    xlabel('x','FontSize',14)
    grid on
hold off

% LTE convergence

figure(2)
lw = 2;
loglog(h_List1,LTE_List,'-.','LineWidth',lw)
hold on
    loglog(h_List,LTE_List1,'-.','LineWidth',lw)
    leg = legend('O(k+h^2)','O(k^2+h^4)');
    set(leg,'FontSize',14);
    title('Convergence','FontSize',14);
    ylabel('log(LTE)','FontSize',14)
    xlabel('log(h)','FontSize',14)
    grid on
hold off

%% Other scheme
%mu = kappa*k/h^2 > 1/6 -->
k = 10*h^2/kappa; 
theta = 1/2 + h^2/(12*kappa*k);
mu = kappa*k/h^2;

%Solving the system
steps = 1000;
Un = solver(m,k,h,theta,x,steps);
Uf = realValue(m,x,steps*k);
U0 = realValue(m,x,0);

%Convergence
m_List = 6:5:200;
k_func = @(h,kappa) 0.3*h^2/kappa;
[LTE_List,h_List,k_List] = LTEconvergence(m_List,kappa,theta,k_func);

%Solving the system
m_List = 6:5:200;
[E_List,h_List1,k_List1] = GlobalError(m_List,kappa,theta);


%%%%% Plotting
%Solution after time
figure(1)
lw = 2;
plot(U0,'-.','LineWidth',lw)
hold on
    plot(Un,'LineWidth',lw)
    plot(Uf,'-.','LineWidth',lw)
    leg = legend('U_0','U_f estimated','U_f real');
    set(leg,'FontSize',14);
    ste = num2str(steps);
    kstring = num2str(k);
    title(['Solution for n = ' ste ' and k = ' kstring],'FontSize',14);
    ylabel('u(x)','FontSize',14)
    xlabel('x','FontSize',14)
    grid on
hold off

% LTE convergence

figure(2)
lw = 2;
loglog(h_List,LTE_List,'-.','LineWidth',lw)
hold on
    title('LTE convergence ','FontSize',14);
    leg = legend('O(k^2+h^2)');
    set(leg,'FontSize',14);
    ylabel('log(LTE)','FontSize',14)
    xlabel('log(h)','FontSize',14)
    grid on
hold off

% Global error Convergence

figure(3)
lw = 2;
loglog(h_List,E_List,'-.','LineWidth',lw)
hold on
    title('Global error after 1 second ','FontSize',14);
    leg = legend('Global error');
    set(leg,'FontSize',14);
    ylabel('log(E)','FontSize',14)
    xlabel('log(h)','FontSize',14)
    grid on
hold off