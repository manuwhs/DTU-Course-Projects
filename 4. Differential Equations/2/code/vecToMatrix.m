function U = vecToMatrix(m,U0)
    U = zeros(m,m);
    for i = 1:m
        U(i,:) = U0((i-1)*m+1:m*i);
    end
end
