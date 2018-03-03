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

h = 0.001;
m = ceil(-1 + 2/h);
x = linspace(-1,1,m); 

k = 0.00005;
tspan = [0 1];
n = ceil((tspan(2)-tspan(1))/k);

%model paramters
epsilon = 0.9/pi;
epsilons = [0.01/pi, 0.1,0.5];
n_graph = 1;
figure()
for ie = 1:n_graph
    epsilon = epsilons(ie);
    %system matrices

    % U0 = sysfunc_3_t0(x, epsilon);
    % U1 = sysfunc_3(x, 1, epsilon);
    % 
    % plot (x, U0);
    % hold on
    % plot (x, U1);
    % 
    % figure();
    % 
    % tmax = 0.1; % Number of seconds

    time_process = 1.6037/pi; % Number of seconds to simulate
    Nsim = ceil(n * time_process); % Number of iterations of the algo
    t = linspace(0,time_process,Nsim);

    %% Set the solution data
    Un = zeros(m,Nsim);
    Un(:,1) = sysfunc_3(x, 0 ,epsilon)';   % t = 0 conditions

    Un(1,:) = sysfunc_3(-1,0, epsilon);   % Boundary conditions
    Un(m,:) = sysfunc_3(1 ,t,epsilon);
    
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

    %% Plotting of the solution
%     subplot(1,n_graph,ie)
%     [T ,X ] = meshgrid(t ,x);
%     mesh( X, T, Un);
%     title(strjoin({'epsilon = ', num2str(epsilon)}))
%     ylabel('time')
%     xlabel('space')
%     zlabel('U(x,t)')
end

 ultimo = Un(:,Nsim);
 Derivatives = diff(ultimo)/h;
 Derivatives = (-ultimo + [ultimo(3:end)' 0 0 ]')/(2*h);
 derivative = Derivatives(floor(m/2));
 