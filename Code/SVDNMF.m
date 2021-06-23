% SVD-NMF initialization

function [W,H] = SVDNMF(A,k);

% This function implements the SVD-NMF algorithm described in [1] for
% initialization of Nonnegative Matrix Factorization Algorithms.
%
% [W,H] = SVDNMF(A,k,flag);
%
% INPUT
% ------------
%
% A    : the input nonnegative m x n matrix A
% k    : the rank of the computed factors W,H
%
% OUTPUT
% -------------
%   
% W   : nonnegative m x k matrix
% H   : nonnegative k x n matrix
%
% 
% References:
% 
% [1] Qiao, H. (2015). New SVD based initialization strategy for 
% non-negative matrix factorization. Pattern Recognition Letters, 63, 71-77.
%
%--------------------------------------------------------------------------


%----------------------check the input matrix------------------------------
if numel(find(A<0)) > 0
    error('The input matrix contains negative elements !')
end
%--------------------------------------------------------------------------

%size of the input matrix
[m,n] = size(A);

% truncated SVD rank-k to the input matrix A. 
[U,S,V] = svds(A,k);
W = abs(U*S); 
H = abs(V)'; 