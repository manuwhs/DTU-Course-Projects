 function [tnList,ynList,hList,rList,nfun] = ImplicitEulersAdaptiveStep(func,...
     Jacob,tspan,N,Y0,abstol,reltol,controller)
%
% This function solves a general first-order Initial Value Problem
% of the form
%                dot_y = f(y,t),  y(tstart) = tbegin
%
% using Implicit euler in n steps (adaptive step size).
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

ynList = [];
tnList = [];
rList = [];
hList = [];
tnList(1) = tbegin;
ynList(:,1)=Y0;

tol = 10e-5; % Newtons method tolerance

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
    
    f = feval(func,tnList(k),ynList(:,k));
    nfun = nfun+1;
    
    acceptedStep = 0;
    it = 0;
    while ~acceptedStep
        it = it+1;
        %y(n+1) in one step
        y_guess = ynList(:,k)+h*f; %Forwards euler
        [y,calls] = newtonODEadaptive(func,Jacob,y_guess,tol,tnList(k)+h,ynList(:,k),h);
        nfun = nfun+calls+1;
        
        %y(n+1) in two step
        hm = h/2;
        ym_guess = ynList(:,k)+hm*f;
        [ym,calls] = newtonODEadaptive(func,Jacob,ym_guess,...
            tol,tnList(k)+hm,ynList(:,k),hm);
        nfun = nfun + calls +1;
        f = feval(func,tnList(k)+hm,ym);
        y_hat_guess = ynList(:,k)+hm*f;
        [y_hat,calls] = newtonODEadaptive(func,Jacob,y_hat_guess,...
            tol,tnList(k)+h,ym,hm);
        nfun = nfun +calls+1;
        
        e = abs(y_hat - y);
        r = max(e./max(abstol,abs(y_hat).*reltol));
        if r<=1.0 || it > 1000
            rList(k) = r;
            hList(k) = h;
            acceptedStep = 1;
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