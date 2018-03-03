function [tnList,ynList,hList,rList,nfun] = RK4AdaptiveStep(func,tspan,N,Y0,abstol,reltol)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                dot_y = f(y,t),  y(tstart) = tbegin
%
% using Rungge Kitta in n steps (adaptive step size).
%
% INPUT:
%    func  : a function handle to function f(u,t)
%    tspan : a 1x2 array of the form [tstart tend]
%    N     : paramter for calculating first step size
%    y0   : initialvalue(s)
%    abstol : absolute tolerance
%    reltol : relative tolerance
% OUTPUT:
%    tnList  : time 
%    ynList : solution
%    hList     : each step size 
%    rList   : estimated error each step
%    nfun : number of function evaluations
%
sizeTspan = size(tspan);
Ninit = sizeTspan(2);  % We can be given only the end time, then the begining is 0

%% Initialization
if (Ninit == 1)
    tbegin = 0;
    tend = tspan;
elseif(Ninit == 2) % [tbegin tend]
    tbegin = tspan(1,1);
    tend = tspan(1,2);
end
h = (tend - tbegin)/N;

ynList = [];
tnList = [];
rList = [];
hList = [];
tnList(1) = tbegin;
ynList(:,1) = Y0;

%%error estimation and control paramters
epstol = 0.8;
facmin = 0.1;
facmax = 5;
kpow = 0.2;

%% Loop
k=1;
nfun = 0;
while tnList(k) < tend
    %In order to compute until tend
    if (tnList(k)+h>tend)
        h = tend-tnList(k);
    end
    f = feval(func,tnList(k),ynList(:,k));
    nfun = nfun+1;
    acceptedStep = 0;
    while ~acceptedStep
        nfun = nfun+10;
        %y(n+1) in one step
        y=RungeKuttaStep(func,tnList(k),ynList(:,k),h,f);
        
        %y(n+1) in two step
        hm = h/2;
        ym=RungeKuttaStep(func,tnList(k),ynList(:,k),hm,f);
        fm = feval(func,tnList(k)+hm,ym);
        y_hat=RungeKuttaStep(func,tnList(k)+hm,ym,hm,fm);
        
        e = abs(y_hat - y);
        r = max(e./max(abstol,abs(y_hat).*reltol));
        if r<=1.0
            acceptedStep = 1;
            hList(k) = h;
            rList(k) = r;
            ynList(:,k+1) = y_hat;
            tnList(k+1) = tnList(k)+h;
        end
        %hnew = h*sqrt(wanted_error/current_error)
        h = max(facmin, min((epstol/r)^kpow,facmax))*h;
    end
    k = k+1;
end