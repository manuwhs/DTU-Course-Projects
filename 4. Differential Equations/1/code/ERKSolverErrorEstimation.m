function [Tout,Xout,Eout] = ...
        ERKSolverErrorEstimation(fun,tspan,x0,h,solver,varargin)
% ERKSOLVERERRORESTIMATION  Fixed step size ERK solver with error est.
%
%                           Solves ODE systems in the form dx/dt = f(t,x)
%                           with x(t0) = x0. 
%
% Syntax:
% [Tout,Xout,Eout]=ERKSolverErrorEstimation(fun,tspan,x0,h,solver,varargin)

% Solver Parameters
s  = solver.stages;     % Number of stages in ERK method
AT = solver.AT;         % Transpose of A-matrix in Butcher tableau
b  = solver.b;          % b-vector in Butcher tableau
c  = solver.c;          % c-vector in Butcher tableau
d  = solver.d;

% Parameters related to constant step size
hAT = h*AT;
hb  = h*b;
hc  = h*c;
hd  = h*d;

% Size parameters
x  = x0;                % Initial state
t  = tspan(1);          % Initial time
tf = tspan(end);        % Final time
N = int64((tf-t)/h);           % Number of steps
nx = length(x0);        % System size (dim(x))

% Allocate memory
T  = zeros(1,s);        % Stage T
X  = zeros(nx,s);       % Stage X
F  = zeros(nx,s);       % Stage F

Tout = zeros(N+1,1);    % Time for output
Xout = zeros(N+1,nx);   % States for output
Eout = zeros(N+1,nx);   % Errors for output

% Algorithm starts here
Tout(1) = t;        
Xout(1,:) = x';
for n=1:N
    % Stage 1
    T(1)   = t;
    X(:,1) = x;
    F(:,1) = fun(T(1),X(:,1),varargin{:});
    
    % Stage 2,3,...,s
    T(2:s) = t + hc(2:s);
    for i=2:s
        X(:,i) = x + F(:,1:i-1)*hAT(1:i-1,i);
        F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
    end

    % Next step
    t = t + h;
    x = x + F*hb;
    e = F*hd;
    
    % Save output
    Tout(n+1) = t;
    Xout(n+1,:) = x';
    Eout(n+1,:) = e';
end