function [yn1_n1,nfun] = newtonODEadaptive(func,Jacob,y0,tol,tn1,yn,h)
%Newton root finding algorithm
%   Fibd the roots of F(yn1) = yn1 - h*f(yn1) - yn 
%   where yn1 is the unknown variable
% INPUT:
%    func  : a function handle to function f(y,t)
%    Jacob : jacobian of f(y,t)
%    y0   :  intial guess for the root i.e yn1_0
%    param : parameters to be passed to func
%    tol : tolerance for the solution
%    yn,tn1,h; are constant paramters of the fucntion F(yn1)
% OUTPUT:
%    yn1_n1  : estimated root

nfun = 0;
yn1_n = y0;
maxIt = 100;
err = 1000;
it = 0;
I = eye(length(y0));

while((it < maxIt) && (norm(err,'inf') > tol))
    it = it +1;
    fprime = I-h*feval(Jacob,tn1,yn1_n);
    f = yn1_n-h*feval(func,tn1,yn1_n)-yn;
    yn1_n1 = yn1_n - fprime\f;
    f = yn1_n1-h*feval(func,tn1,yn1_n1)-yn;
    nfun = nfun+2;
    err = f;
    yn1_n = yn1_n1;
end
end

