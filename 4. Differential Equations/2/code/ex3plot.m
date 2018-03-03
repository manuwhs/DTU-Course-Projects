hstring = num2str(h);
lw = 2;
figure(1)
hold on
    surf(X,Y,Uhat)
    surf(X,Y,Usol)
    title(['Estimated Solution for h = ' hstring],'FontSize',14);
    xlabel('x');
    ylabel('y');
hold off
figure(2)
hold on
    title(['Real Solution for h = ' hstring],'FontSize',14);
    surf(X,Y,Usol)
hold off

if(calculateError==1)
figure(3)
loglog(hList,Elist,'-.','LineWidth',lw)
hold on
    loglog(hList,Elist9_deferred,'--','LineWidth',lw)
    title('Global Error convergence','FontSize',14); 
    leg = legend('9-point stencil','9- points Deferred solution');
    set(leg,'FontSize',14);
    ylabel('log(E)','FontSize',12,'FontWeight','bold')
    xlabel('log(h)','FontSize',12,'FontWeight','bold')
    grid on
hold off
end

% figure(3)
% hold on
%     plot(x,sol(x,y(20)))
%     plot(x,Umatrix(20,:))
%     legend('real','estimation')
% hold off
% 
% figure(4)
% hold on
%     plot(x,sol(x(1),y))
%     plot(x,Umatrix(:,1))
%     legend('real','estimation')
% hold off