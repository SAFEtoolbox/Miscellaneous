% This script shows some examples of how to handle
% discrete uniform distributions when using the sampling functions
% 'AAT_sampling.m' and 'OAT_sampling.m' of the SAFE Toolbox
%
% This script prepared by F. Pianosi               
% University of Bristol, 2017                                           
% mail to: francesca.pianosi@bristol.ac.uk                                     
                                                                                                                      
                                                                               
% Consider a case of M=2 input factors (x1,x2), 
% x1 must be sampled from a continous uniform distribution in [0,2]
% x2 must be sampled from a discrete uniform distribution in [1,5]

% The SAFE code to generate a Latin Hypercube of N=10 samples of x1,x2 is:
SampStrategy = 'lhs' ;
N = 10 ; 
M = 2  ;
DistrFun = {'unif','unid'} ;
DistrPar = { [ 0 2 ], 5 } ;
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N) ;

% Notice that samples of x2 (i.e. the second column of X):
X(:,2)
% vary between 1 and 5 even if in DistrPar{2} we only specified 
% the upper bound of the range (5). 
% This is because AAT_sampling uses the Matlab/Octave function
% 'unidinv' function, which assumes by default that the lower bound be 1 
% (hence it does not need specifying).

% What if we want to sample from a range with lower bound different from 1?
% We suggest to still use the function AAT_sampling and shift the
% results after sampling. 

% For example, say we want to sample x2 from a discrete uniform 
% distribution in [0,5], then we should use AAT_sampling to sample 
% from [1,6] and then subtract 1, i.e.:
DistrPar = { [ 0 2 ], 6 } ;
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N) ;
X(:,2) = X(:,2) - 1 ;
X(:,2) % Check the results

% Other examples:
% Code to sample x2 from a discrete uniform in [10,15]:
DistrPar = { [ 0 2 ], 6 } ;
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N) ;
X(:,2) = X(:,2) + 9 ;
X(:,2)
% Code to sample x2 from a discrete uniform in [-3,0]:
DistrPar = { [ 0 2 ], 4 } ;
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N) ;
X(:,2) = X(:,2) - 4 ;
X(:,2)

% The same applies for OAT_sampling. 
% For example, let's implement the last case above (x2 in [-3,0]):
r = 10 ;
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,'trajectory');
X(:,2) = X(:,2) - 4 ;
X(:,2)

