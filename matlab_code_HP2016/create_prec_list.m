function Xprec = create_prec_list(prec,rM,Nr)
%
% Xprec = create_prec_list(prec,rM,Nr)
%
% prec = time series to be perturbed                    - vector (Nsteps,1)
%   rM = allowed ratio of variability (e.g. rM=0.4 means +/-40%)   - scalar
%   Nr = number of perturbed time series to be generated           - scalar
%
% Xprec = ensemble of perturbed time series            - matrix (Nsteps,Nr)

Nsteps = length(prec);
Xprec = nan(Nsteps,Nr) ;
for i=1:Nr
    Xprec(:,i) = perturbe_rainfall(prec,1-rM,1+rM) ;
end
%figure
%plot(Xprec,'Color',[126 126 126]/256)
%hold on; plot(prec, 'k')

function [rain_perturbed,rMult] = perturbe_rainfall(rain,rMult_min,rMult_max)
%
% [rain_perturbed,rMult] = perturbe_rainfall(rain,rMult_min,rMult_max)

% Find timesteps where a rainfall event (storm) begins and ends:
rain_start = find( diff([0; rain>0; 0])==1 );
rain_end   = find( diff([0; rain>0; 0])==-1 ) ;
if rain_end(end)>length(rain); rain_end(end)=length(rain); end

% Perturb rainfall time series using a random multiplier for each storm:
rain_perturbed = rain ;
rMult = nan(size(rain_start));
for i=1:length(rain_start)
    rMult(i) = rMult_min + rand*(rMult_max-rMult_min) ;
    rain_perturbed(rain_start(i):rain_end(i)) = rain(rain_start(i):rain_end(i))*rMult(i) ;
end