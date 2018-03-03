function [Jacob] = JacobDepPrey( t, Y )
%DEPPREY Summary of this function goes here
%   Detailed explanation goes here

%     a = 1;
%     b = 1;
    x = Y(1);
    y = Y(2);

    % Syntax: xdot = PreyPredator(t,x,a,b)
    Jacob = zeros(2);
    Jacob(1,1) = 1-y;
    Jacob(2,1) = y;
    Jacob(2,2) = -1+x;
    Jacob(1,2) = -x;
end
