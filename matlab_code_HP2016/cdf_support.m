function [YF, Y_min,Y_max] = cdf_support(YU,YC,threshold,N)

[M,n]       = size(YC) ;
[NU,Nsteps] = size(YU) ;

% find minimum and maximum value of 'y' across all samples
% (i.e. the unconditional sample YU and all conditional subsamples in 'YC')
Y_min = min(YU) ;
Y_max = max(YU) ;
for i=1:M; 
    for k=1:n ;
        Y_min = min([Y_min;YC{i,k}]);
        Y_max = max([Y_max;YC{i,k}]);
    end
end

YF=cell(1,Nsteps) ;
for t=1:Nsteps
    DYt = (threshold(t)-Y_min(t))/(N-1)  ; 
    YF{t} = [ Y_min(t):DYt:threshold(t) Y_max(t) ] ;
end


