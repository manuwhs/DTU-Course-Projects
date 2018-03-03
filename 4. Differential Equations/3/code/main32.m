%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Exercise 3: Hyperbolic equations   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subscripts
addpath('functions');

clc;
clear all;
close all;

%%


%model paramters
epsilon = 0.01;
epsilons = [0.01/pi, 0.1,0.5];
n_graph = 1;

Errors = {};
Error_abs = [];
nplots = 0;
n_errors_seen = 1;

hs = [0.02, 0.01, 0.0075, 0.005, 0.002,0.001];
nshow = 5;
for ih = 1:nshow
    % time and space grids

    h = hs(ih);
    m = ceil(-1 + 2/h);
    x = linspace(-1,1,m); 

    k =  epsilon * h^2;
    
    tspan = [0 1];
    n = ceil((tspan(2)-tspan(1))/k);
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

    time_process = 2*k; % Number of seconds to simulate
    Nsim = ceil(n * time_process); % Number of iterations of the algo

    
    t = linspace(0,time_process,Nsim);

    n_errors_seen = 1;
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
    exact_solution = sysfunc_3(x,tn, epsilon);
    
    % Compute error
    Error = exact_solution' - Un(:,end);
    Errors{n_errors_seen} = Error;
    Error_abs = [Error_abs norm(Error)];
    n_errors_seen = n_errors_seen + 1;
end
figure()

plot(log(hs(1:nshow)), log(Error_abs), 'color', rand(1, 3), 'DisplayName', 'k = epsilon * h *h', 'LineWidth', 3)

xlabel('log(h)' ,'FontSize', 50)
ylabel('log(LTE)', 'FontSize', 50)
legend('show')
hold off
grid on
 ultimo = Un(:,Nsim);
 Derivatives = diff(ultimo)/h;
 Derivatives = (-ultimo + [ultimo(3:end)' 0 0 ]')/(2*h);
 derivative = Derivatives(floor(m/2));
 