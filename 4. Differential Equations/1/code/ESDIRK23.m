function [Tout,Xout,Gout,Eout,info,stats] = ESDIRK23(fun,jac,tspan,x0,h0,absTol,relTol,Method,varargin)

%% ESDIRK23 Parameters 
%=========================================================================
% Runge-Kutta method parameters
Gout = 0;

gamma = 1-1/sqrt(2);
a31 = (1-gamma)/2;
AT = [0 gamma a31;0 gamma a31;0 0 gamma];
c  = [0; 2*gamma; 1];
b  = AT(:,3);
bhat = [    (6*gamma-1)/(12*gamma); ...
    1/(12*gamma*(1-2*gamma)); ...
    (1-3*gamma)/(3*(1-2*gamma))    ];
d  = b-bhat;
p  = 2;
phat = 3;
s = 3;

% error and convergence controller
epsilon = 0.8;
tau = 0.1*epsilon; %0.005*epsilon;
itermax = 20;

%========================================================================

t0 = tspan(1);
tf = tspan(2);
info = struct(...
            'Method',    Method,  ... % carsten
            'nStage',    s,       ... % carsten
            'absTol',    'dummy',  ... % carsten
            'relTol',    'dummy',  ... % carsten
            'iterMax',   itermax, ... % carsten
            'tspan',     tspan,   ... % carsten
            'nFun',      0, ...
            'nJac',      0, ...
            'nLU',       0, ...
            'nBack',     0, ...
            'nStep',     0, ...
            'nAccept',   0, ...
            'nFail',     0, ...
            'nDiverge',  0, ...
            'nSlowConv', 0);


        
%% Main ESDIRK Integrator
%========================================================================
nx = size(x0,1);
F = zeros(nx,s);
t = t0;
x = x0;
h = h0;
I = eye(nx);

F(: ,1) = feval (fun ,t,x, varargin {:}) ;

 info . nFun = info . nFun +1;
 dfdx = feval (jac ,t,x, varargin {:}) ;
 info . nJac = info . nJac +1;
 if (t+h)>tf
 h = tf -t;
 end
 hgamma = h* gamma ;
 dRdx = I - hgamma * dfdx ;
 [L,U, pivot ] = lu(dRdx ,'vector');
 info . nLU = info . nLU +1;

 iter = zeros (1,s);
 

% Output. Initial size of the output, it will increase if more step than
% chunk are needed.
chunk = 100;
Tout = zeros(chunk,1);
Xout = zeros(chunk,nx);
Eout = zeros(chunk,nx); 

Tout(1,1) = t;
Xout(1,:) = x.';
Eout(1,:) = zeros(1,nx); 

% While we have not reached the end of the tspan
while t<tf
    

    info.nStep = info.nStep+1;
    %=====================================================================
    % A step in the ESDIRK method
    i=1;   % Variables with the step number in each iteration of the algo.
    diverging = false;
    SlowConvergence = false; % carsten
    alpha = 0.0;
    Converged = true;
    
    % While we have not reached the last step of the method and the
    % Newton's method has converged.
    
    while (i<s) && Converged 
        
   
        % Stage i=2,...,s of the ESDIRK Method
        i=i+1;
        phi = x + F(: ,1:i -1) *(h*AT (1:i -1,i));

        % Initial guess for the state

        dt = c(i)*h;
        X = x + dt*F(: ,1);
        T = t+dt;
            
        F(:,i) = feval (fun ,T,X, varargin {:}) ;
        info.nFun = info.nFun+1;
        R = X - hgamma*F(:,i) - phi;
%        rNewton = norm(R./(absTol + abs(G).*relTol),2)/sqrt(nx);
        rNewton = norm(R./(absTol + abs(X).*relTol),inf);
        Converged = (rNewton < tau);
        %iter(i) = 0; % original, if uncomment then comment line 154: iter(:) = 0;
        
        
        % Newton Iterations !!!
        % Apply Newton's method iterations until it converges, for every
        % implicit step of the ESDIRK model
        while ~Converged && ~diverging && ~SlowConvergence%iter(i)<itermax
            iter(i) = iter(i)+1;
            dX = U\(L\(R(pivot,1)));
            info.nBack = info.nBack+1;
            X = X - dX;
            rNewtonOld = rNewton;
            F(:,i) = feval(fun,T,X,varargin{:});
            info.nFun = info.nFun+1;
            R = X - hgamma*F(:,i) - phi;
%            rNewton = norm(R./(absTol + abs(G).*relTol),2)/sqrt(nx);
            rNewton = norm(R./(absTol + abs(X).*relTol),inf);
            alpha = max(alpha,rNewton/rNewtonOld);
            Converged = (rNewton < tau);
            diverging = (alpha >= 1);
            SlowConvergence = (iter(i) >= itermax); % carsten
            %SlowConvergence = (alpha >= 0.5); % carsten
            %if (iter(i) >= itermax), i, iter(i), Converged, diverging, pause, end % carsten
        end
        % diverging will have the first diverging state !! 
        % It it diverges, it will get out of the previous loop.
        %diverging = (alpha >= 1); % original, if uncomment then comment line 142: diverging = (alpha >= 1)*i;
        diverging = (alpha >= 1)*i; % carsten, recording which stage is diverging
    end
    %if diverging, i, iter, pause, end
    nstep = info.nStep;
    stats.t(nstep) = t;
    stats.h(nstep) = h;
    stats.r(nstep) = NaN;
    stats.iter(nstep,:) = iter;
    stats.Converged(nstep) = Converged;
    stats.Diverged(nstep)  = diverging;
    stats.AcceptStep(nstep) = false;
    stats.SlowConv(nstep)  = SlowConvergence*i; % carsten, recording which stage is converging to slow (reaching maximum no. of iterations)
    iter(:) = 0; % carsten
    
    
    %=====================================================================
    % Error and Convergence Controller
  
    % If all the implicit staged converged !!
    
    if Converged
        % Error estimation
        e = F*(h*d);
        
%        r = norm(e./(absTol + abs(G).*relTol),2)/sqrt(nx);
        r = norm(e./(absTol + abs(X).*relTol),inf);
        r = max(r,eps);
        stats.r(nstep) = r;
       
        % Next Step
        t = T;
        x = X;
        F(:,1) = F(:,s);            

        % Just taking care that we do not go further than tspan
        h = max(1e-8,h);
        if (t+h) > tf
            h = tf-t;
        end
         % Jacobian Update Strategy
         dfdx = feval (jac ,t,x, varargin {:}) ;
         info . nJac = info . nJac +1;
         hgamma = h* gamma ;
         dRdx = I - hgamma * dfdx ;
         [L,U, pivot ] = lu(dRdx ,'vector');
         info . nLU = info . nLU +1;
         
    % If we were not able to converge in one of the steps...
    else % not converged
        info.nFail=info.nFail+1;
        CurrentStepAccept = false;

        hgamma = h*gamma;
        dRdx = dgdx - hgamma*dfdx;
        [L,U,pivot] = lu(dRdx,'vector');
        info.nLU = info.nLU+1;
        hLU = h;
    end
    
    %=====================================================================
    % Storage of variables for output
       info . nAccept = info . nAccept + 1;
       nAccept = info.nAccept;
       if nAccept > length(Tout);
           Tout = [Tout; zeros(chunk,1)];
           Xout = [Xout; zeros(chunk,nx)];
           Eout = [Eout; zeros(chunk,nx)];
           
       end
       Tout(nAccept,1) = t;
       Xout(nAccept,:) = x.';
       Eout(nAccept,:) = e;
        
end
info.nSlowConv = length(find(stats.SlowConv)); % carsten
nAccept = info.nAccept;
Tout = Tout(1:nAccept,1);
Xout = Xout(1:nAccept,:);
Eout = Eout(1:nAccept,:);
