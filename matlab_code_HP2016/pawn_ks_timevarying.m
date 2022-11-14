function [KS_max, KS, counter, YF, FU ] = pawn_ks_timevarying(YU,YC,YF,threshold)


[M,n] = size(YC);
[NU,Nsteps] = size(YU) ;

%% Determine support for CDF estimation

% % For each time step, pre-allocate the support to the unconditional output
% % sample:
% YF=cell(1,Nsteps) ; for t=1:Nsteps; YF{t} = YU(:,t) ; end 
% % Then, for each input and conditioning value, add any other lower or
% % higher value from the conditional samples:
% for i=1:M; for k=1:n ; for t=1:Nsteps
%             YCik=YC{i,k}; YF{t} = unique([YF{t};min(YCik(:,t));max(YCik(:,t))]); 
% end; end; end
% % figure; hold on; for t=1:Nsteps; plot(t,YF{t},'ok'); end

%% Estimate unconditional CDF (for each time step)

FU=cell(1,Nsteps) ;
for t=1:Nsteps
    FU{t} = empiricalcdf(YU(:,t),YF{t});
end

%% Compute KS statistic
FC     = cell(1,Nsteps) ;
KS     = cell(M,1)      ;
counter= cell(M,1)      ;
KS_max =  nan(M,Nsteps) ;
for i=1:M
    KSi = nan(n,Nsteps) ;
    counteri = nan(n,Nsteps) ;
    for k=1:n
        YCik = YC{i,k} ;
        for t=1:Nsteps
            FC{t} = empiricalcdf(YCik(:,t),YF{t});
            FUt = FU{t} ;
            FCt = FC{t} ;
            
            % compute KS over subrange:
            idx = YF{t} <= threshold(t) ;
            if sum(idx)==0; idx = true(size(YF{t})) ; fprintf('i=%d,k=%d,t=%d:no output samples below threshold, KS computed over avialable range of output\n',i,k,t); end
            KSi(k,t) = max(abs(FUt(idx)-FCt(idx))) ;
            counteri(k,t) = sum(idx) ;
            
            % Compute KS over entire range:
            %KSi(k,t) = max(abs(FU{t}-FC{t}))     ;
            
        end       
    end
    counter{i}  = counteri ;
    KS{i}       = KSi      ;
    KS_max(i,:) = max(KSi) ;
    fprintf('Input %d completed\n',i)
end
