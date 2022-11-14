function [ X_EET_radial, Y_EET_radial] = fromVBSAtoEET(XB,YB,XC,YC)
%
% [ X_EET_radial, Y_EET_radial] = fromVBSAtoEET(XB,YB,XC,YC)
%
% This functions takes the input/output samples (XB,YB,XC,YC) generated
% to compute the VBSA first-order and total-order indices 
% and re-arrange them into the matrices (X_EET_radial,Y_EET_radial) that
% can be used to compute EET indices (with 'radial' design).
%
% For more details and references on the sampling/resampling strategy 
% to obtain (XB,YB,XC,YC) and approximation of VBSA indices, 
% see help of functions: 
% VBSA_resampling, VBSA_indices
%
% For more details and references on the sampling strategy (including 
% radial design) and definition of EET indices,
% see help of functions:
% OAT_sampling, EET_indices
%
% Input:
%    XB = set of input samples                             - matrix (N,M)
%    XC = set of input samples from resampling             - matrix (N*M,M)
%    YB = set of output samples                            - matrix (N,P)
%    YC = set of output samples from resampling            - matrix (N*M,P)
%
% Output:
% X_EET_radial = set of input samples ready for EET    - matrix (N*(M+1),M)
% Y_EET_radial = set of output samples ready for EET   - matrix (N*(M+1),P)

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%


[NX,M]=size(XC) ;
[NY,P]=size(YC) ;
if NX~=NY; error('Input ''XC'' and ''YC'' must have the same number of rows.'); end
N = NX/M ;
if abs(N-round(N)); error('The number of rows in input ''XC'' must be a multiple of the number of columns.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute X_EET_radial and Y_EET_radial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = N ; % the number of EEs (and blocks in matrix 'X_EET_radial'), which
% we usually denote by 'r', is equal to the base sample size (N) used for
% VBSA
X_EET_radial = nan(r*(M+1),M);
Y_EET_radial = nan(r*(M+1),P);
for j=1:N 
    idx = (j-1)*(M+1)+1; % row at which current block (i.e. block j) starts
    % (for example, j=1 -> idx=1; j=2 -> idx=M+1+1; etc.)
    X_EET_radial(idx,:) = XB(j,:) ; % first component of the block comes from XB
    Y_EET_radial(idx,:) = YB(j,:) ;
    for i=1:M % all subsequent components of the block will come from XC
        X_EET_radial(idx+i,:) = XC((i-1)*N+j,:) ;
        Y_EET_radial(idx+i,:) = YC((i-1)*N+j,:) ;
    end
end

% A bit more explanation of how the code works:
%
% Matrix AB that would be produced by OAT_sampling can be obtained
% by rearranging rows of XA and XB produced by AAT_sampling, i.e.:
%
% AB = reshape([XB(:) XA(:)]',2*N, M); % AB = (2*N,M)
% % = [ XB(1,:)
% %     XA(1,:)
% %     XB(2,:)
% %     XA(2,:)
% %      ...
% %     XB(N,:)
% %     XA(N,:) ]
%
% Matrix X that would be produced from AB by OAT_sampling 
% is instead the matrix X_EET_radial produced by this function.


