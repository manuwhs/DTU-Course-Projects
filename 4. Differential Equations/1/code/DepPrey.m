function [ gradY] = DepPrey( t, Y, args)
%DEPPREY Summary of this function goes here
%   Detailed explanation goes here

    a = 1;
    b = 1;

    % Syntax: xdot = PreyPredator(t,x,a,b)
    gradY = zeros(2,1);
    gradY(1) = a*(1-Y(2))*Y(1);
    gradY(2) = -b*(1-Y(1))*Y(2);
end

