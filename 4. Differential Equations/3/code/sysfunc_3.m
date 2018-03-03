function [ out ] = sysfunc_3( x,t,epsilon)
     out = -tanh((x - 0.5 -t )/(2*epsilon)) +1;
     % out = -sin(pi * x );
end

