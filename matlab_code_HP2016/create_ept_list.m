function Xept = create_ept_list(ept,eM,Ne)
%
% Xept = create_ept_list(ept,eM,Ne)
%
% ept = time series to be perturbed                     - vector (Nsteps,1)
%  eM = allowed ratio of variability (e.g. eM=0.4 means +/-40%)    - scalar
%  Ne = number of perturbed time series to be generated            - scalar
% Xept = ensemble of perturbed time series             - matrix (Nsteps,Ne)

Nsteps = length(ept);
Xept = nan(Nsteps,Ne) ;
for i=1:Ne
    Xept(:,i) = ept*( 1-eM+2*eM*rand ) ;
end

%figure
%plot(Xept,'Color',[126 126 126]/256)
%hold on; plot(ept, 'k')
