function Xflow = create_flow_list(flow,eQ,alfa,Nf)
%
% Xflow = create_flow_list(flow,eQ,alfa,Nf)
%
% flow = time series to be perturbed                    - vector (Nsteps,1)
%   eQ = allowed ratio of variability (in probabilistic terms, 
%        e.g. eQ=0.4 means +/-40% with probability 0.99)           - scalar
% alfa = lag-1 autocorrelation of errors                           - scalar
%   Nf = number of perturbed time series to be generated           - scalar
%
% Xflow = ensemble of perturbed time series            - matrix (Nsteps,Nf)


beta = eQ/3*sqrt(1-alfa^2) ;
Nsteps = length(flow);
Xflow = nan(Nsteps,Nf) ;
for i=1:Nf
    Xflow(:,i) = perturbe_flow(flow,alfa,beta) ;
end
%figure
%subplot(211); plot(Xflow,'Color',[126 126 126]/256); hold on; plot(flow, 'k','LineWidth',2)
%subplot(212); plot(Xflow-repmat(flow,1,Nf),'Color',[126 126 126]/256); hold on; plot(zeros(size(flow)), 'k')
%
%tmp = ( Xflow./repmat(flow,1,Nf) < 1+eQ )+(Xflow./repmat(flow,1,Nf) > 1-eQ)==2 ;
%figure; plot( sum(tmp)/Nsteps,'o')

function [flow_perturbed,pert] = perturbe_flow(flow,alfa,beta)
%
% [flow_perturbed,pert] = perturbe_flow(flow,alfa,beta)

sigma_err = beta*flow ; 

pert = nan(size(flow));
pert(1)=0;
for t=1:length(flow)-1
    pert(t+1)=alfa*pert(t)+ randn*sigma_err(t) ;
    if flow(t+1)+pert(t+1)<0 ; pert(t+1)=0; end
end
flow_perturbed = flow + pert ;