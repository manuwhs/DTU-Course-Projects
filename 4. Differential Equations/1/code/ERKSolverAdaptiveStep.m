function [tnList,ynList,Eout] = ...
        ERKSolverAdaptiveStep(fun,tspan,Y0,h,solver,abstol,reltol,varargin)
% ERKSOLVERERRORESTIMATION  Fixed step size ERK solver with error est.
%
%                           Solves ODE systems in the form dx/dt = f(t,x)
%                           with x(t0) = x0. 
%
% Syntax:
% [Tout,Xout,Eout]=ERKSolverErrorEstimation(fun,tspan,x0,h,solver,varargin)

% Solver Parameters
s  = solver.stages;     % Number of stages in ERK method
AT = solver.AT;         % Transpose of A-matrix in Butcher tableau
b  = solver.b;          % b-vector in Butcher tableau
c  = solver.c;          % c-vector in Butcher tableau
d  = solver.d;

% Parameters related to constant step size
hAT = h*AT;
hb  = h*b;
hc  = h*c;
hd  = h*d;

% Size parameters
x  = Y0;                % Initial state
t  = tspan(1);          % Initial time
tend = tspan(end);        % Final time
N = int64((tend-t)/h);           % Number of steps
nx = length(Y0);        % System size (dim(x))

% Allocate memory
T  = [];        % Stage T
X  = [];       % Stage X
F  = [];       % Stage F

tnList = [];    % Time for output
ynList = [];   % States for output
Eout = [];   % Errors for output

%%error estimation and control paramters
epstol = 0.8;
facmin = 0.1;
facmax = 5;
kpow = 1/6;

% Algorithm starts here
tnList(1) = t;        
ynList(:,1) = x';
k = 1;
while tnList(k)<tend
    
    if (tnList(k)+h>tend)
            h = tend-tnList(k);
    end
    
    acceptedStep = 0;
    while ~acceptedStep
        
        hAT = h*AT;
        hb  = h*b;
        hc  = h*c;
        hd  = h*d;   
        
        % Stage 1
        T(1)   = t;
        X(:,1) = x;
        F(:,1) = fun(T(1),X(:,1),varargin{:});

        % Stage 2,3,...,s
        T(2:s) = t + hc(2:s);
        for i=2:s
            X(:,i) = x + F(:,1:i-1)*hAT(1:i-1,i);
            F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
        end

        % Next step
        t = t + h;
        x = x + F*hb;
        e = F*hd;
        r = max(e./max(abstol,abs(x).*reltol));
        if k == 2
            k
        end
        if r<=1.0
            acceptedStep = 1;
%             hList(k) = h;
%             rList(k) = r;
            ynList(:,k+1) = x;
            tnList(k+1) = tnList(k)+h;
%             Eout(n+1,:) = e';
          
        end
        %hnew = h*sqrt(wanted_error/current_error)
        h = max(facmin, min((epstol/r)^kpow,facmax))*h;
        
    end
    k = k+1;
    tnList(k)
    k
end
% while tnList(k) < tend
%     
%     if (tnList(k)+h>tend)
%             h = tend-tnList(k);
%     end
% 
%     acceptedStep = 0;
%     while ~acceptedStep
%         %In order to compute until tend
%         
%         hAT = h*AT;
%         hb  = h*b;
%         hc  = h*c;
%         hd  = h*d;   
%         % Stage 1
%         T(1)   = t;
%         X(:,1) = x;
%         F(:,1) = fun(T(1),X(:,1),varargin{:});
%         % Stage 2,3,...,s
%         T(2:s) = t + hc(2:s);
%         for i=2:s
%             X(:,i) = x + F(:,1:i-1)*hAT(1:i-1,i);
%             F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
%         end% Next step
%         x_full = x + F*hb;
% 
%         %%%%%%%%%%% 2 steps %%%%%%%%%%%%%
%         hAT = (h/2)*AT;
%         hb  = (h/2)*b;
%         hc  = (h/2)*c;
%         hd  = (h/2)*d;   
%         % Stage 1
%         T(1)   = t;
%         X(:,1) = x;
%         F(:,1) = fun(T(1),X(:,1),varargin{:});
%         % Stage 2,3,...,s
%         T(2:s) = t + hc(2:s);
%         for i=2:s
%             X(:,i) = x + F(:,1:i-1)*hAT(1:i-1,i);
%             F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
%         end% Next step
%         t = t + (h/2);
%         x_half = x + F*hb;
% 
%         %%%%%%%%%%%%%  
%         % Stage 1
%         T(1)   = t;
%         X(:,1) = x_half;
%         F(:,1) = fun(T(1),X(:,1),varargin{:});
%         % Stage 2,3,...,s
%         T(2:s) = t + hc(2:s);
%         for i=2:s
%             X(:,i) = x_half + F(:,1:i-1)*hAT(1:i-1,i);
%             F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
%         end% Next step
%         t = t + (h/2);
%         x_hat = x_half + F*hb;
%         
%         e = abs(x_hat - x_full);
%         r = max(e./max(abstol,abs(x_hat).*reltol));
%         if r<=1.0
%             acceptedStep = 1;
% %             hList(k) = h;
% %             rList(k) = r;
%             ynList(:,k+1) = x_hat;
%             tnList(k+1) = tnList(k)+h;
% %             Eout(n+1,:) = e';
%         end
%         %hnew = h*sqrt(wanted_error/current_error)
%         h = max(facmin, min(sqrt(epstol/r),facmax))*h;
%         
%     end
%     %%%%%%%%%%
%     k=k+1;
% end