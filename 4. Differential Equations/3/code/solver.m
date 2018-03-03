function [Un] = solver(m,k,h,theta,x,steps)

%system matrices
[Al,Ae] = sysMatrices(m,k,h,theta);
U0 = realValue(m,x,0);

%Solving the system
Un = U0;
for i = 1:steps
    tn = k*(i-1);
    tn1 = k*i;
    gn = boundaries(m,x(1),x(end),tn);
    gn1 = boundaries(m,x(1),x(end),tn1);
    %calculate next step
    Un1 = Al\(Ae*Un + (gn1-gn));
    
    %update
    Un = Un1;
end

end

