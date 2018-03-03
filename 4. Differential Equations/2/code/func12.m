function f = func12(x,u,epsilon)
%DEPPREY Summary of this function goes here
%   Detailed explanation goes here
    
    % Syntax: xdot = PreyPredator(t,x,a,b)
    f = zeros(2,1);
    f(1) = u(2);
    f(2) = (u(1)-u(1)*u(2))/epsilon;
end

