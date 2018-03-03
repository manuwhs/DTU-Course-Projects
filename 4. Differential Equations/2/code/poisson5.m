function A = poisson5(m)
    e = ones(m,1);
    S = spdiags([e -2*e e], [-1 0 1], m, m);
    I = speye(m);
    A = kron(I,S)+kron(S,I);
    A=(m+1)^2*A;
end