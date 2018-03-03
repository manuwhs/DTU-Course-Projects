function [E_List,h_List,k_List] = GlobalError(m_List,kappa,theta)

E_List = zeros(1,length(m_List));
h_List = zeros(1,length(m_List));
k_List = zeros(1,length(m_List));
time = 1;

for i = 1:length(m_List)
    h = 2/(m_List(i)-1);
    k = 0.3*h^2/kappa;        
    x = linspace(-1,1,m_List(i));
    steps = round(time/k);
    Un = solver(m_List(i),k,h,theta,x,steps);
    Uf = realValue(m_List(i),x,steps*k);
    E_List(i) = norm(Un-Uf);
end
end