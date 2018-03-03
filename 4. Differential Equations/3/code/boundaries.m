function g = boundaries(m,x0,xend,t)
%%% Calculate boundary values

%system constants
kappa =0.1;
Q = 3;
alpha = [1,4,16];
ai = 1;
bi = 0;

%bounraies 
gr = 0;
gl = 0;
for i=1:Q
    gl = gl + (ai*cos(alpha(i)*x0)+bi*sin(alpha(i)*x0))* ...
        exp(-kappa*alpha(i)^2*t);
    gr = gr + (ai*cos(alpha(i)*xend)+bi*sin(alpha(i)*xend))* ...
        exp(-kappa*alpha(i)^2*t);
end

g = zeros(m,1);
g(1) = gl;
g(end) = gr;

end