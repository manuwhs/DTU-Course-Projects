%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Exercise 2.3: Hyperbolic equations   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subscripts
addpath('functions');

clc;
clear all;
close all;

%%
% time and space grids

h = 0.01;
m =  + 2/h;
x = linspace(-1,1,m); 

k = h*1.6;
tspan = [0 1];
n = (tspan(2)-tspan(1))/k;

%model paramters
a = 0.5;
c = a;
f = 0;

%system matrices

U0 = sin(2*pi*x);

%Solving the system
Un = U0;

A1 = 1 - (c*k/h)*(1 - exp(sqrt(-1)* 2* pi * h));
A1 = A1^1;
A1_mod = abs(A1);
A1_phase = angle(A1);
phase_ideal = 2 * pi * k * a;

A1_f = A1^(40*n);
A1_fmod = abs(A1_f );
A1_fphase = angle(A1_f );
phase_ideal_f = rem(2 * pi * k * a * (40*n), 2\pi);
total_delat = A1_fphase - phase_ideal_f;

figure(1)
% pl = plot(x,Un, 'color',[0 0 1],'LineStyle','-', 'LineWidth', 3);
% hold on 
% exact_solution = sin(2*pi*(x -c*0));
% pl = plot(x,exact_solution, 'color',[1 0 0],'LineStyle','--', 'LineWidth', 3);
% hold on 

legend('Estimated','Real')
xlabel('x (space)' ,'FontSize', 50)
ylabel('u(x,t)', 'FontSize', 50)
Errors = {};
Error_abs = [];
nplots = 1;
n_errors_seen = 1;

for i = 1:80*n
    tn = k*(i-1);
    tn1 = k*i;
    % i
    %calculate next step
    Un1 = Un;
    N_elem = length(Un1);
    for j = 2:(N_elem-1)
        Un1(j) = Un(j) - (c*k/h)*(Un(j) - Un(j-1));
    end
    Un1(1) = Un(1) - (c*k/h)*(Un(1) - Un(N_elem));
    Un1(N_elem) = Un(N_elem) - (c*k/h)*(Un(N_elem) - Un(N_elem-1));
    
    %Un1(N_elem) = Un1(1)
    
    %update
    Un = Un1;
    if i >= nplots*(n)*80
%         subplot(3,1,nplots + 1);

   %     Ux = A1_fmod * sin(2*pi*x + A1_fphase);
        hold on 
        plot(x,Un1, 'color',[nplots/8 nplots/8 1] ,'LineStyle','-', 'LineWidth', 3)
        
        hold on
        exact_solution = sin(2*pi*(x -c*tn1));
        plot(x,exact_solution,'color',[1 nplots/8 nplots/8],'LineStyle','--', 'LineWidth', 3)
        nplots = nplots + 1;

        % Compute error
        Error = abs(exact_solution - Un1);
        Errors{n_errors_seen} = Error;
        Error_abs = [Error_abs norm(Error)];
        n_errors_seen = n_errors_seen + 1;
    end
end

plot ([-1 1],[0,0])

A1 = 1 - (c*k/h)*(1 - exp(sqrt(-1)* 2* pi * h));
A1 = A1^1;
A1_mod = abs(A1);
A1_phase = angle(A1);
phase_ideal = 2 * pi * k * a;

A1_f = A1^(80*n);
A1_fmod = abs(A1_f );
A1_fphase = angle(A1_f );
phase_ideal_f = rem(2 * pi * k * a * (40*n), 2\pi);
total_delat = A1_fphase - phase_ideal_f;


delay = A1_phase - phase_ideal;
A2 = 1 - (c*k/h)*(1 - exp(-sqrt(-1)* 2* pi * h));
A2 = A2^1;
A2_mod = abs(A2);
A2_phase = angle(A2);

den = exp(sqrt(-1)* 2* pi * h) - exp(- sqrt(-1)* 2* pi * h);
num = A1 * exp(sqrt(-1)* 2* pi * h) - A2 * exp(- sqrt(-1)* 2* pi * h);
num2 = exp(sqrt(-1)* 2* pi * (h - a*k)) -  exp(- sqrt(-1)* 2* pi * (h - a*k));

Total_gain = num / den; 
Total_gain_mod = norm(Total_gain);
Total_gain_angle = angle(Total_gain);


figure()
plot(sin(2*pi*x - 0.5))
hold on
plot(sin(2*pi*x))
