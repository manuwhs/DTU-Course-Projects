function F = form_RHS5(m,f,bound,x,y)
%Return right-hand side of a 2D multidimensional linear system
%AU = F where A m^2xm^2 matrix and F = m^2x1 vector
%f(i,j) = f(x(i),y(j)) = F(m*i+j)
    h = 1/(m+1);
    x0 = x(1)-h; xend = x(end)+h;
    y0 = y(1)-h; yend = y(end)+h;
    F = zeros(m*m,1);
    for j = 1:m
        for i =1:m
            F((j-1)*m+i) = feval(f,x(i),y(j));
        end
        %Dirilecht BC in x
        F((j-1)*m+1) =F((j-1)*m+1) -feval(bound,x0,y(j))/h^2;
        F(j*m) = F(j*m) -feval(bound,xend,y(j))/h^2;
    end
    %Dirilecht BC in y
    for i = 1:m
        F(i) = F(i) - feval(bound,x(i),y0)/h^2;
    end
    for i = 1:m
        F(m*(m-1)+i) = F(m*(m-1)+i) - feval(bound,x(i),yend)/h^2;
    end
end