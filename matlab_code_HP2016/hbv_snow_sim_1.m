function [TSsim,TSobs] = hbv_snow_sim_1(X,Xprec,temp,Xept,Xflow,Xparsnow,Xparsoil,Xparroute,Case)

N      = size(X,1)     ;
Nsteps = size(Xprec,1) ;
TSsim  = nan(N,Nsteps) ; 
TSobs  = nan(N,Nsteps) ;
for j=1:N
    prec  = Xprec(:,X(j,1)) ;
    ept   = Xept(:,X(j,2))  ;    
    P     = snow_routine(temp,prec,Xparsnow(X(j,3),:)) ;
    Q_sim = hbv_sim(P,ept,[Xparsoil(X(j,4),:),Xparroute(X(j,5),:)],Case);
    TSsim(j,:) = Q_sim'           ;
    TSobs(j,:) = Xflow(:,X(j,6))' ;    
end