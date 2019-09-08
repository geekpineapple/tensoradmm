
function [x,history] = tensor_cpl_admm(A,y,rho,alpha,sX,maxItr,varargin)

%set Default
ABSTOL        =    1e-5                                ;
RELTOL        =    1e-5                                ;
t_start       =    tic                                 ;

%% ADMM solver
[~,n]         =    size(A)                             ;
x             =    zeros(n,1)                          ;
z             =    zeros(n,1)                          ;
u             =    zeros(n,1)                          ;
zhat          =    zeros(n,1)                          ;
uhat          =    zeros(n,1)                          ;
q             =    sparse(y)                           ;

for k         =     1 : maxItr   
    zold      =     z                                  ;
    uold      =     u                                  ;
    one       =     ones(n,1)                          ;
    one       =     diag(sparse(double(one)))          ;   
    x         =     (A'*A+rho*one)\(rho*zhat - uhat + A'*q)          ;    
    z         =     reshape((x + 1/rho*uhat),sX);
    z         =     tvdenoise(z,rho,3)                 ; 
    z         =     tvdenoise(z,0.003,3)               ;
    z         =     tvdenoise(z,0.001,3)               ;
    z         =     tvdenoise(z,0.0007,3)              ;
    z         =     z(:)                               ;
    u         =     uhat + rho*(x - z)                 ;
    r_norm    =     norm(x-z)                          ;
    s_norm    =     norm(-rho*(z - zold))              ;
    %% apply acceleration.
    Ek = -1;
    if k > 1
       Ek = max(history.r_norm(k-1),...
               history.s_norm(k-1))-max(r_norm,s_norm) ; 
    end
   
    if Ek > 0
        alpha_old = alpha                              ;
        alpha =  (1+sqrt(1+4*alpha^2))/2               ;
        zhat  =  z + (alpha_old-1)/alpha*(z - zold)    ;
        uhat  =  u + (alpha_old-1)/alpha*(u - uold)    ;
        
    else
        alpha =  1                                     ;
        uhat  =  u                                     ;
        zhat  =  z                                     ; 
    end   
    %% diagnostics, reporting, termination checks
    history.r_norm(k)   =  r_norm                      ;
    history.s_norm(k)   =  s_norm                      ; 
    history.eps_pri(k)  = sqrt(n)*ABSTOL + ...
                          RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k) = sqrt(n)*ABSTOL + ...
                          RELTOL*norm(rho*u)           ;
    fprintf('%3d\t%10.4f\t%10.4f\t%10.6f\t%10.6f\t\n', k, ...
           history.r_norm(k), history.eps_pri(k), ...
           history.s_norm(k), history.eps_dual(k))      ;

    if (history.r_norm(k) < history.eps_pri(k) && ...
        history.s_norm(k) < history.eps_dual(k))
        break;
    end
end
toc(t_start);
end