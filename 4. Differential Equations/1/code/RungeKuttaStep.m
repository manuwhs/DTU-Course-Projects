function [y1] = RungeKuttaStep(func,t,y,h,f)
%SINGLE STEP RUNGE-KUTTA
s1 = f;  % Explicit Euleer
s2 = feval(func,t + h/2 ,y + (h/2) * s1);  % Midpoint method
s3 = feval(func,t + h/2 ,y + (h/2) * s2);  % Recursive of Midpoint method
s4 = feval(func,t + h,y + s3*h); % Slope in the predictive point as if used Euler with s3 as slope

y1 = y+ (h/6)*(s1 + 2*s2 + 2*s3 + s4);

end

