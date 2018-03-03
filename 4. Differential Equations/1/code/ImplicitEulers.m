 function [tnList,ynList] = ImplicitEulers(func,Jacob,tspan,N,Y0)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                Ydot = F(y,t),  y(tbegin) = Y0
%
% using Implicit Euler Method with adaptive time step 
%
% INPUT:
%    func   : a function handle to function F(Y,t)
%    Jacob  : a function handle to function dF(Y,t)/dY
%    tspan  : a 1x2 array of the form [tbegin tend]
%    N      : total number of steps in tspan
%    Y0     : initialvalue(s)s
%

sizeY = size(Y0);
Ndim = sizeY(1);
sizeTspan = size(tspan);
Ninit = sizeTspan(2);

%% Initialization
if (Ninit == 1)
    tbegin = 0;
    tend = tspan;
elseif(Ninit == 2) % [tbegin tend]
    tbegin = tspan(1,1);
    tend = tspan(1,2);
end
h = (tend - tbegin)/N;

ynList = zeros(Ndim,N+1);
tnList = zeros(1,N+1);
tol = 10e-5; % Newtons method tolerance
ynList(:,1)=Y0;

%% Loop
for k = 1:N
    f = feval(func,tnList(k),ynList(:,k));
    y_guess = ynList(:,k)+h*f; %Forwards euler
    ynList(:,k+1) = newtonODE(func,Jacob,y_guess,tol,tnList(k)+h,ynList(:,k),h);
    tnList(:,k+1) = tnList(:,k)+h;
end