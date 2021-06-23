% SVD-based initialization for NMF referred to as 
% Nonnegative SVD with low-rank correction (NNSVDLRC); see 
% Improved SVD-based Initialization for Nonnegative Matrix Factorization 
% using Low-Rank Correction, Syed Muhammad Atif, Sameer Qazi, Nicolas
% Gillis, July 2018. 

function [W,H,Y,Z,e] = NNSVDLRC(X,r,delta,maxiter);

% *********
%   Input
% *********
% X      : input nonnegative m-by-n matrix.
% r      : rank of the sought NMF decomposition.
% --- Optional ---
% delta  : stoppong parameter for the low-rank correction (default: 0.05).
% maxiter: maximum number of iteration for the low-rank correction 
%           (default: 20).
%
% **********
%   Output
% **********
% (W,H)  : m-by-r and r-by-n nonnegative matrices s.t. WH approximate X.
% (Y,Z)  : rank-p truncated SVD of X so that YZ approximates X, 
%          where p = ceil(r/2+1).
% e      : relative error ||X-WH||_F/||X||_F of the initialization (W,H)
%          throughout the low-rank correction updates. 

if nargin <= 2
    delta = 0.05;
end
if nargin <= 3
    maxiter = 20;
end
p = floor(r/2+1);
if issparse(X)
    [u,s,v] = svds(X,p); 
else
    [u,s,v] = mySVD(X,p); 
end
%     [u,s,v] = svds(X,p); 

Y = u*sqrt(s); 
Z = sqrt(s)*v'; 
% Best rank-one approximation
W(:,1) = abs(Y(:,1)); 
H(1,:) = abs(Z(1,:)); 
% Next (r-1) rank-one factors
i = 2; j = 2; 
while i <= r
    if mod(i,2) == 0
        W(:,i) = max(Y(:,j),0); 
        H(i,:) = max(Z(j,:),0); 
    else
        W(:,i) = max(-Y(:,j),0); 
        H(i,:) = max(-Z(j,:),0);
        j = j+1;
    end
    i = i+1; 
end
% Scale (W,H): this is important for HALS. 
WtYZ = (W'*Y)*Z; 
WtW = W'*W; 
HHt = H*H'; 
scaling = sum(sum(WtYZ.*H))/sum(sum((WtW).*(HHt))); 
H = H*sqrt(scaling); 
W = W*sqrt(scaling); 
WtYZ = WtYZ*sqrt(scaling); 
WtW = WtW*scaling; 
HHt = HHt*scaling; 
% relative error
nX = sqrt( sum ( sum( (Y'*Y).*(Z*Z') ) ) ); 
e(1) = sqrt( nX^2 - 2*sum(sum( WtYZ.*H ) ) + sum(sum( WtW.*(HHt) )) ) / nX; 
% Improve using LRA-based accelerated HALS
% until relative error improvment is below delta of the initial error. 
k = 1; 
while (k == 1 || e(k-1)-e(k) > delta*e(1)) && k <= maxiter 
    W = LRAnnlsHALSupdt(Z',Y',H',W')'; 
    [H,WtW,WtX] = LRAnnlsHALSupdt(Y,Z,W,H); 
    e(k+1) = sqrt( nX^2 - 2*sum(sum( WtX.*H ) ) + sum(sum( WtW.*(H*H') )) ) / nX; 
    k = k+1; 
end