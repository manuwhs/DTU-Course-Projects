%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Exercise 2.2: Hyperbolic equations   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subscripts
addpath('functions');

clc;
clear all;
close all;

%%
% time and space grids

figure(1)

m_list = [50, 70, 80, 100, 120,200]; % space 

n_list = [140, 160, 180, 200,300]; % time

for in = 1:length(n_list) % time

    Error_abs = [];
    k_s = [];
    h_s = [];

   for im = 1:length(m_list)  % space
        n = n_list(1 + length(n_list) - in);
        tspan = [0 1];
        k = (tspan(2)-tspan(1))/n;
        k_s = [ k_s k];

        m = m_list(1 + length(m_list) -im);
        h = 2/(m-1);
        h_s = [ h_s h];

        x = linspace(-1,1,m); 

        %model paramters
        a = 0.5;
        c = a;
        f = 0;

        %system matrices

        U0 = sin(2*pi*x);

        %Solving the system
        Un = U0;

        nplots = 1;
        n_errors_seen = 1;
        for i = 1:200   % 200
            tn = k*(i-1);
            tn1 = k*i;
            gn = boundaries(m,x(1),x(end),tn);
            gn1 = boundaries(m,x(1),x(end),tn1);

            %calculate next step
            Un1 = Un;
            N_elem = length(Un1);
            for j = 2:(N_elem-1)
                Un1(j) = Un(j) - (c*k/h)*(Un(j) - Un(j-1));
            end
            Un1(1) = Un(1) - (c*k/h)*(Un(1) - Un(N_elem-1));
            Un1(N_elem) = Un(N_elem) - (c*k/h)*(Un(N_elem) - Un(N_elem-1));

            %update
            Un = Un1;

        end
        exact_solution = sin(2*pi*(x -c*tn1));
        
        Error = abs(exact_solution - Un1);

        Error_abs = [Error_abs max(Error)];
    end

    hold on 
    
    %% Check stability 
    
    caca = sum(k > 2*h_s);
    if (caca == 0)  
        plot(h_s, Error_abs, 'color', rand(1, 3), 'DisplayName', strjoin({'k = ', num2str(k)}), 'LineWidth', 3)
    end
end
xlabel('m (space-delta)' ,'FontSize', 50)
ylabel('max(|error|)', 'FontSize', 50)
legend('show')
hold off
    

