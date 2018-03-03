function [Al,Ae] = sysMatrices(m,k,h,theta)
%%% Calculate system matrices

if(theta<0 || theta>1)
    disp('Invalid theta');
    return
end


%model paramter
kappa = 0.1;
mu = k*kappa/h^2;

%A0 matrix construction
mid = -2.*ones(m,1);
sides = ones(m,1);
A0 = spdiags([sides mid sides],[-1,0,1],m,m);
A0(1,1) = 0; A0(1,2) = 0; A0(end,end) = 0; A0(end,end-1) = 0;

I = speye(m);
Al = (I-mu*theta*A0);
Ae = (I+mu*(1-theta)*A0); 


end