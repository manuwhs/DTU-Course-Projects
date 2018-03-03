function [tnList,ynList,hList,rList,nfun] = ExplicitEulersAdaptiveStep(func,...
    tspan,N,Y0,abstol,reltol,controller)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                dot_y = f(y,t),  y(tstart) = tbegin
%
% using Explicit Euler in n steps (adaptive step size).
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

sizeY = size(Y0);
Ndim = sizeY(1);
sizeTspan = size(tspan);
Ninit = sizeTspan(2);  % We can be given only the end time, then the begining is 0

ynList = [];
tnList = [];
rList = [];
hList = [];

%% Initialization
if (Ninit == 1)
    tbegin = 0;
    tend = tspan;
elseif(Ninit == 2) % [tbegin tend]
    tbegin = tspan(1,1);
    tend = tspan(1,2);
end
h = (tend - tbegin)/N;

tnList(1) = tbegin;
ynList(:,1) = Y0;
%%error estimation and control paramters
epstol = 0.8;
facmin = 0.1;
facmax = 5;

if controller == 'PI'
    a = 1/6;
    b = 1/6;
    c = 1/2;
else
    a = 1/2;
    b = 0;
    c = 0;
end

%% Loop
k=1;
nfun = 0;
while tnList(k) < tend
    
    %In order to compute until tend
    if (tnList(k)+h>tend)
        h = tend-tnList(k);
    end
    nfun = nfun+1;
    f = feval(func,tnList(k),ynList(:,k));
    
    acceptedStep = 0;
    while ~acceptedStep
        nfun = nfun+1;
        y = ynList(:,k)+h*f;
        
        hm = h/2;
        ym = ynList(:,k)+hm*f;
        f = feval(func,tnList(k)+hm,ym);
        y_hat = ym + hm*f;
        
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
        if k == 1
            con = ((epstol/r)^a);
        else
            con = ((epstol/r)^a)*((epstol/rList(k-1))^b)*((h/hList(k-1))^-c);
        end
        h = max(facmin, min(con,facmax))*h;
    end
    k = k+1;
end