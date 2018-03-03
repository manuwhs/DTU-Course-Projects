function A = poisson9(m)
%This function return the 9-points stencil shceme for the discrete
% Poisson equation
    e = ones(m,1);
    S = spdiags([e e e], [-1 0 1], m, m);
    C = spdiags([e 4*e e], [-1 0 1], m, m);
    B = spdiags([3*e -24*e 3*e], [-1 0 1], m, m);
    I = speye(m);
    A = kron(S,C)+kron(I,B);
    A = (m+1)^2*A./6;
end