function U0 = realValue(m,x,t)
%%% Calculate real u(x,t) values

%system constants
kappa =0.1;
Q = 3;
alpha = [1,4,16];
ai = 1;
bi = 0;

U0 = zeros(m,1);

for j=1:m
    for i=1:Q
        U0(j) = U0(j) + (ai*cos(alpha(i)*x(j))+bi*sin(alpha(i)*x(j)))* ...
            exp(-kappa*alpha(i)^2*t);
    end
end

end