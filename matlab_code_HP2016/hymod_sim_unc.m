function [RMSE,TS] = hymod_sim_unc(X,Prec,Ept,ParSoil,ParRoute,Flow)
%
% 
% Prec

N      = size(X,1)    ;
Nsteps = size(Prec,1) ;
TS     = nan(N,Nsteps); 
TSobs  = nan(N,Nsteps);
for j=1:N
    prec  = Prec(:,X(j,1)) ;
    ept   = Ept(:,X(j,2))  ;
    param = [ ParSoil(X(j,3),:) ParRoute(X(j,4),:) ]  ;
    Q_sim = hymod_sim(param,prec,ept);
    Q_obs = Flow(:,X(j,5)) ;
    TS(j,:)   = Q_sim' ;
    TSobs(j,:)= Q_obs' ;    
end

% RMSEt = tv_objfun(TS(:,warmup+1:end),TSobs(:,warmup+1:end),'RMSE',W) ;
RMSE  = sqrt(mean((TS-TSobs).^2,2));