function [EList,hList] = globalError5(mList,func,sol,span)
%GLOBAL_ERROR this function return the global error for each mList
%INPUTS
% mList: list of number of points
% func: function handle of RHS  AU = F
% sol:  function handle for the ture U solution
% x:    space discretization vector
% y:    space discretization vector
% span: grid size
%OUTPUS:
% Elist: global error for each mList
% hlist: step size for each mList

    n = length(mList);
    EList = zeros(n,1);
    hList = zeros(n,1);
    for i = 1:n
       m = mList(i);
       h =  abs(span(1)-span(2))/(m+1);
       x = linspace(span(1)+h,span(2)-h,m);
       y = linspace(span(1)+h,span(2)-h,m);
       %Create scheme
       A = poisson5(m);
       F = form_RHS5(m,func,sol,x,y);
       %Get solution
       Uvec = A\F;
       Uhat = vecToMatrix(m,Uvec);
       
       %true solution
       [X,Y]=meshgrid(x,y); 
       Usol = sol(X,Y);
       
       %Global error
       EList(i) = h*norm(Usol-Uhat);
       hList(i) = h;
    end
end

