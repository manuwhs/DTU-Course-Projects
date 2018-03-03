function [ Uk ] = newtonODE(functions,U0,h,tol,epsilon)
%Newton root finding algorithm
%   Find the roots of non-linear set of functions G(u)
% INPUT:
%    functions  : function handle for G and J
%    h :    step size
%    U0   :  vector of length m with intial guess for every u(x)
%    tol :  tolerance for the solution
% OUTPUT:
%    Uk  : estimated root
    Uk = U0;
    [G,J] = functions(Uk,h,epsilon);
    maxIt = 100;
    it = 0;
    while((it < maxIt) && (norm(G,'inf') > tol))
        it = it +1;      
        Uk(2:end-1) = Uk(2:end-1) - J\G;
        [G,J] = functions(Uk,h,epsilon);
    end
end

