function rmse = RMSE(y_sim,y_obs)
%
% Computes the Root Mean Squared Error
%
% rmse = RMSE(Y_sim,Y_obs)
% 
% Y_sim = time series of modelled variable     - matrix (N,T)
%         (N>1 different time series can be evaluated at once)
% Y_obs = time series of observed variable     - matrix (N,T)
%         (N==1: the same time series of observations is used
%          to compare with each time series in Y_sim)
%
% rmse  = vector of RMSE coefficients          - vector (N,1)

[N,T] = size(y_sim) ;
[M,D] = size(y_obs) ;

if T~=D
    error('input y_sim and y_obs must have the same number of columns')
end
if M==1
    y_obs = repmat(y_obs,N,1) ;
elseif M~=N
    error('input y_sim and y_obs must have the same number of rows')
end

Err  = y_sim - y_obs ;
rmse = sqrt(mean(Err.^2,2)) ;