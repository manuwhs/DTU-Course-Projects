function [LTE_List,h_List,k_List] = LTEconvergence(m_List,kappa,theta,kfunc)

LTE_List = zeros(1,length(m_List));
h_List = zeros(1,length(m_List));
k_List = zeros(1,length(m_List));
steps = 1;

for i = 1:length(m_List)
    h = 2/(m_List(i)-1);
    k = kfunc(h,kappa);       
    x = linspace(-1,1,m_List(i));
    Un = solver(m_List(i),k,h,theta,x,steps);
    Uf = realValue(m_List(i),x,steps*k);
    LTE_List(i) = norm(Un-Uf);
    h_List(i) = h;
    k_List(i) = k;
end
end

