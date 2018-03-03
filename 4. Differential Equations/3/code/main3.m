%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Exercise 3: Hyperbolic equations   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subscripts
addpath('functions');

clc;
clear all;
close all;

%%
% time and space grids

h = 0.01;
m = ceil(-1 + 2/h);
x = linspace(-1,1,m); 

k = 0.2*h^2/epsilon;

%model paramters
epsilon = 0.5;

time_process = 1.5; % Number of seconds to simulate
Nsim = round(time_process/k);
t = zeros(1,Nsim);
for i = 1:Nsim
    t(i) = k*i;
end

% Set the solution data
Un = zeros(m,Nsim);
% func1 = @(x) -sin(pi*x);
% Un(1,:) = zeros(1,Nsim);
% Un(m,:) = zeros(1,Nsim);
% Un(:,1) = func1(x);
func = @(x,t,epsilon) -tanh(x+0.5-t/(2*epsilon))+1;
Un(1,:) = func(-1,t,epsilon);
Un(m,:) = func(1 ,t,epsilon);
Un(:,1) = func(x,0,epsilon);
Uf = func(x,t(end),epsilon);

for i = 2:Nsim
    tn = k*(i-1);
    tn1 = k*i;
    %calculate next step
    for j = 2:(m-1)
        Un(j,i) = (k*epsilon/h^2 )*(Un( j-1,i-1) - ...  %% 1
            2*Un(j, i - 1) + Un(j + 1, i-1)) + Un(j, i - 1) -...  %% 2
             k/h*Un(j,i-1) * (Un(j,i-1)-Un(j-1,i-1));   %% 3
    end
end

figure(1)
lw = 2;
plot(Un(:,1),'--','LineWidth',lw)
hold on
    plot(Un(:,end),'-.','LineWidth',lw)
    plot(Uf,'--','LineWidth',lw)
    legend('Initial','Final Solution','Final Solution Real')
    title('Solution of the grid','FontSize',14);
    ylabel('U(x)','FontSize',14)
    xlabel('x','FontSize',14)
    grid on
hold off