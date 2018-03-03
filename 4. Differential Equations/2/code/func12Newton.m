function f = func12Newton(x,u,epsilon)
%DEPPREY Summary of this function goes here
%   Detailed explanation goes here
    
    % Syntax: xdot = PreyPredator(t,x,a,b)
    f = zeros(4,1);
    f(1) = u(2);
    f(2) = (u(1)-u(1)*u(2))/epsilon;
    f(3) = u(4);
    f(4) = (u(1)-1)*u(3)/epsilon - u(1)*u(4)/epsilon;
end

