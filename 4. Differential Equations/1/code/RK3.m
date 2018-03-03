function [Tout,Xout,Eout] = RK3(func,tspan,Y0,h, varargin)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                u? = f(u,t),  u(tstart) = eta
%
% using Euler?s Method in n steps (constant step size).
%
% INPUT:
%    func  : a function handle to function f(u,t)
%    tspan : a 1x2 array of the form [tstart tend]
%    N     : total number of steps in tspan
%    y0   : initialvalue(s)
%    param : parameters to be passed to func
%

%% Define the parameters of the model
A = [0 0 0;
    1/2 0 0;
    -1 2 0];
b = [1/6; 2/3; 1/6];
c = [0; 1; 1/2];
d = [1/12; -1/6; 1/12];
s = 3;

AT = A.';
% Parameters related to constant step size
hAT = h*AT;
hb  = h*b;
hc  = h*c;
hd  = h*d;

% Size parameters
x  = Y0;                % Initial state
t  = tspan(1);          % Initial time
tf = tspan(end);        % Final time
N = int64((tf-t)/h);           % Number of steps
nx = length(Y0);        % System size (dim(x))

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
    F(:,1) = func(T(1),X(:,1),varargin{:});
    
    % Stage 2,3,...,s
    T(2:s) = t + hc(2:s);
    for i=2:s
        X(:,i) = x + F(:,1:i-1)*hAT(1:i-1,i);
        F(:,i) = feval(func,T(i),X(:,i),varargin{:});
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