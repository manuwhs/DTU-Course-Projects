function F = form_RHS9(m,f,bound,lapla,x,y,referred)
%Return right-hand side of a 2D multidimensional linear system
%AU = F where A m^2xm^2 matrix and F = m^2x1 vector
%f(i,j) = f(x(i),y(j)) = F(m*i+j)
    h = 1/(m+1);
    %Dirichlet BC
    x0 = x(1)-h; xend = x(end)+h;
    y0 = y(1)-h; yend = y(end)+h;
    u0j = bound(x0,[y0 y yend]);        %BC vector along u(x(0),y(j))
    uendj = bound(xend,[y0 y yend]);    %BC vector along u(x(end),y(j))
    ui0 = bound([x0 x xend],y0);        %BC vector along u(x(i),y(0))
    uiend = bound([x0 x xend],yend);    %BC vector along u(x(i),y(end))
    
    F = zeros(m*m,1);
    for j = 1:m
        for i =1:m
            F((j-1)*m+i) = feval(f,x(i),y(j));
        end
        %Dirilecht BC in x;
        F((j-1)*m+1) = F((j-1)*m+1) - 4*u0j(j+1)/(6*h^2) - u0j(j+2)/(6*h^2) ...
            - u0j(j)/(6*h^2);
        F(j*m) = F(j*m) - 4*uendj(j+1)/(6*h^2) - uendj(j)/(6*h^2) - ...
            uendj(j+2)/(6*h^2);
    end
    % Dirilecht BC in y
    % boundary conditions along y(0) 
    for i = 1:m
        F(i) = F(i) - 4*ui0(i+1)/(6*h^2) - ui0(i)/(6*h^2) - ui0(i+2)/(6*h^2);
    end
    % corners u(0,0) and u(1,0) have beaing added twice. We substract once
    F(m) = F(m) + ui0(end)/(6*h^2);
    F(1) = F(1) + ui0(1)/(6*h^2);
    
    % boundary conditions along y(end) 
    for i = 1:m
        F(m*(m-1)+i) = F(m*(m-1)+i) -4*uiend(i+1)/(6*h^2) -uiend(i)/(6*h^2) ...
            -uiend(i+2)/(6*h^2);
    end
    % corners u(1,1) and u(0,1) have beaing added twice. We substract once
    F(m*m) = F(m*m) + uiend(end)/(6*h^2);
    F(m*(m-1)+1) = F(m*(m-1)+1) + uiend(1)/(6*h^2);
    
    if(referred==1)
        for j = 1:m
            for i =1:m
                F((j-1)*m+i) = F((j-1)*m+i)+(h^2/12)*feval(lapla,x(i),y(j));
            end
        end   
    end
end