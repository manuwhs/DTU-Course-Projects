function [G,J] = functions11(u,h,epsilon)
%%% Caluclate the multivariate non-linear funcation G(u) and its Jacobian
% J(u)

G = epsilon/h^2*( u(1:end-2) - 2*u(2:end-1) + ...
    u(3:end))+ u(2:end-1).*((u(3:end) - u(1:end-2))/(2*h) - 1);

J = 1/h^2 * (diag(-2*epsilon + 0.5*h*(u(3:end) - u(1:end-2))-h^2) ...
    + diag(epsilon - 0.5*h*u(3:end-1),-1) + ...
    diag(epsilon + 0.5*h*u(2:end-2),1));

end