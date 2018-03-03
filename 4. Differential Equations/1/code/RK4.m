function [tnList,ynList] = RK4(func,tspan,N,Y0)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                dot(Y) = F(Y,t),  Y(tstart) = Y0
%
% using classic Runge Kutta method.
%
% INPUT:
%    func  : a function handle to function f(y,t)
%    tspan : a 1x2 array of the form [tstart tend]
%    N     : total number of steps in tspan
%    Y0   : initialvalue(s)
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
    s1 = feval(func,tnList(k),ynList(:,k));  % Explicit Euleer
    s2 = feval(func,tnList(k) + dt/2 ,ynList(:,k) + (dt/2) * s1);  % Midpoint method
    s3 = feval(func,tnList(k) + dt/2 ,ynList(:,k) + (dt/2) * s2);  % Recursive of Midpoint method
    s4 = feval(func,tnList(k) + dt,ynList(:,k) + s3*dt); % Slope in the predictive point as if used Euler with s3 as slope

    ynList(:,k+1) = ynList(:,k)+ (dt/6)*(s1 + 2*s2 + 2*s3 + s4);
    tnList(k+1) = tnList(k)+dt;
end