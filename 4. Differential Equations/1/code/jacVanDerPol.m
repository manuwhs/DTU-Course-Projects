function [ jac ] = jacVanDerPol( t, Y, args )
%VANDELPOL Summary of this function goes here
%   Detailed explanation goes here
    mu = 100;
    x = Y(1);
    y = Y(2);
    jac = zeros(length(Y));
    % Syntax: xdot = PreyPredator(t,x,a,b)
    jac(1,1) = 0;
    jac(1,2) = 1;
    jac(2,1) = -2*mu*x*y -1;
    jac(2,1) = mu*(1-x^2);
end

