%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Exercise 2: 9-points Laplacian  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add paths
addpath('functions');
%% Parametes
close all
m = 500; 
span = [0 1];
h = abs(span(2)-span(1))/(m+1);
x = linspace(span(1)+h,span(2)-h,m);
y = linspace(span(1)+h,span(2)-h,m);
%RHS and boundary functions

sol31 = @(x,y) sin(4*pi*(x+y))+cos(4*pi*x.*y);
func31 = @(x,y) (4*pi)^2*(-2*sin(4*pi*(x+y))-cos(4*pi*x*y)*(x^2+y^2));
lapla31 = @(x,y) (4*pi)^2*(64*pi^2*sin(4*pi*(x+y))-4*cos(4*pi*x*y)+ ...
    32*pi*x*y*sin(4*pi*x*y)+16*pi^2*cos(4*pi*x*y)*(x^2+y^2)^2);
functions1 = {func31,sol31,lapla31};

sol32 = @(x,y) x.^2 + y.^2;
func32 = @(x,y) 4;
lapla32 = @(x,y) 0;
functions2 = {func32,sol32,lapla32};

sol33 = @(x,y) sin(2*pi*abs(x-y).^(2.5));
func33 = @(x,y) 0;
lapla33 = @(x,y) 0;
functions3 = {func33,sol33,lapla33};

% sol33 = @(x,y) sin(2*pi*(abs(x-y).^(2.5)));
% func33 = @(x,y)  2 * -5*pi*(x^2-2*y*x+y^2)* ...
%     (10*pi*abs(x-y)^(5/2)*sin(2*pi*abs(x-y)^(5/2))...
%     -3*cos(2*pi*abs(x-y)^(5/2)))/(2*abs(x-y)^(3/2));
% lapla33 = @(x,y) 2 * -5*pi*((abs(x-y)*...
%     (80*pi*abs(x-y)^(13/2)+(380*pi*x^2-760*pi*y*x+380*pi*y^2)...
%     *abs(x-y)^(9/2)+(-70*pi*x^4+280*pi*y*x^3-420*pi*y^2*x^2+....
%     280*pi*y^3*x-70*pi*y^4)*abs(x-y)^(5/2))+...
%     (-1000*pi^3*x^10+10000*pi^3*y*x^9- ...
%     45000*pi^3*y^2*x^8+120000*pi^3*y^3*...
%     x^7-210000*pi^3*y^4*x^6+252000*pi^3*y^5*x^5-...
%     210000*pi^3*y^6*x^4+120000*pi^3*y^7*x^3-45000*...
%     pi^3*y^8*x^2+10000*pi^3*y^9*x-1000*pi^3*y^10)*...
%     abs(x-y)^(5/2))*sin(2*pi*abs(x-y)^(5/2))+...
%     ((3*x^4-12*y*x^3+18*y^2*x^2-12*y^3*x+3*y^4)...
%     *abs(x-y)+1800*pi^2*x^10-18000*pi^2*y*x^9+...
%     81000*pi^2*y^2*x^8-216000*pi^2*y^3*x^7+378000*...
%     pi^2*y^4*x^6-453600*pi^2*y^5*x^5+378000*pi^2*y^6*...
%     x^4-216000*pi^2*y^7*x^3+81000*pi^2*y^8*x^2-18000 ...
% *pi^2*y^9*x+1800*pi^2*y^10)*cos(2*pi*abs(x-y)^(5/2)))/(8*abs(x-y)^(13/2));
% functions3 = {func33,sol33,lapla33};

%% 9-points Laplacian on function 1
close all
calculateError = 0;             %calculate and plot Error convergence
plotting = 1;                   %show plots 
corrected = 0;                  %deferred corrections
funcs = functions3;
%Create scheme and RHS
A = poisson9(m);
F = form_RHS9(m,funcs{1},funcs{2},funcs{3},x,y,corrected);
%Get estimated solution
Uvec = A\F;
Uhat = vecToMatrix(m,Uvec);
%Get real solution
[X,Y]=meshgrid(x,y);
Usol = funcs{2}(X,Y);
E = h*norm(Usol-Uhat); %Global error

%Global error convergence
if(calculateError==1)
    mList = [4,5,8,10,15,20,30,50,100,200,500]; 
    [Elist9,hList] = globalError(mList,funcs{1},funcs{2},...
        funcs{3},span,corrected);
end
%Global error convergence withouf deferred correction
if(calculateError==1 && corrected ==1)
    [Elist9_deferred,hList] = globalError(mList,funcs{1}, ...
        funcs{2},funcs{3},span,0);
end
%Plotting
if(plotting==1)
    Elist = Elist9;
    ex3plot
end

%%  5-points stencil
close all
calculateError = 1;             %calculate and plot Error convergence
plotting = 1;                   %show plots
funcs = functions3;

%Create scheme and RHS
A = poisson5(m);
F = form_RHS5(m,funcs{1},funcs{2},x,y);
%Get estimated solution
Uvec = A\F;
Uhat = vecToMatrix(m,Uvec);
%Get real solution
[X,Y]=meshgrid(x,y);
Usol = funcs{2}(X,Y);
E = h*norm(Usol-Uhat); %Global error

%Global error convergence
if(calculateError==1)
    mList = [4,5,8,10,15,20,30,50,100,200,500]; 
    [Elist5,hList] = globalError5(mList,funcs{1},funcs{2},span);
end
%Plotting
if(plotting==1)
    Elist = Elist5;
    ex3plot
end


