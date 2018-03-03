function [tnList,ynList] = ExplicitEulers(func,tspan,N,Y0)
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
sizeY = size(Y0);
Ndim = sizeY(1);
sizeTspan = size(tspan);
Ninit = sizeTspan(2);  % We can be given only the end time, then the begining is 0

ynList = zeros(Ndim,N+1);
tnList = zeros(1,N+1);

%% Initialization
if (Ninit == 1)
    tbegin = 0;
    tend = tspan;
elseif(Ninit == 2) % [tbegin tend]
    tbegin = tspan(1,1);
    tend = tspan(1,2);
end
dt = (tend - tbegin)/N;

tnList(1) = tbegin;
ynList(:,1) = Y0;

%% Loop

for k = 1:N
    f = feval(func,tnList(k),ynList(:,k));
    ynList(:,k+1) = ynList(:,k)+dt*f;
    tnList(k+1) = tnList(k)+dt;
end