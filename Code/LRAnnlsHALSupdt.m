function [V,UtU,UtM] = LRAnnlsHALSupdt(Y,Z,U,V,alphaparam,delta)

% Optimizing the following nonnegative least squares (NNLS) problem  
% 
%    min_{V >= 0} ||M-UV||_F^2, where M = Y*Z, 
% 
% with an exact block-coordinate descent scheme, updating the rows of V
% sequentially using a closed-form solution. 
% 
% Code adapted from 
% Gillis, N., & Glineur, F. (2012). Accelerated multiplicative updates 
% and hierarchical ALS algorithms for nonnegative matrix factorization. 
% Neural computation, 24(4), 1085-1105.
% 
% *********
%   Input
% *********
% (Y,Z) : input m-by-p and p-by-n matrices. 
%         Note that YZ is not necessarily nonnegative.  
% (U,V) : matrix U that defines the NNLS problem, and initialization for V
% 
% --- Optional ---
% alphaparam and delta control the number of iterations (default: 0.5, 0.1,
% respectively). For more information,  we refer the reader to 
% Gillis, N., & Glineur, F. (2012). Accelerated multiplicative updates 
% and hierarchical ALS algorithms for nonnegative matrix factorization. 
% Neural computation, 24(4), 1085-1105. 
%
% **********
%   Output
% **********
% V    : r-by-n matrix that is has lower error ||XY - UV|| than the initial
%        matrix V provided in input. 
% UtU  : = U'*U. 
% UtM  : = U'*M. 

if nargin <= 4
    alphaparam = 0.5; 
end
KX = sum( Y(:) > 0 ) + sum( Z(:) > 0 );
n = size(V,2); 
[m,r] = size(U);
maxiter = floor( 1 + alphaparam*(KX+m*r)/(n*r+n) );
if nargin <= 5
    delta = 0.01; 
end
%% Precomputations 
UtU = U'*U;
UtM = (U'*Y)*Z; 
%% Coordinate descent 
eps0 = 0; cnt = 1; eps = 1; 
while eps >= (delta)^2*eps0 && cnt <= maxiter
    nodelta = 0; if cnt == 1, eit3 = cputime; end
    cputime1 = cputime; 
        for k = 1 : r
            deltaV = max((UtM(k,:)-UtU(k,:)*V)/UtU(k,k),-V(k,:));
            V(k,:) = V(k,:) + deltaV;
            nodelta = nodelta + deltaV*deltaV'; % used to compute norm(V0-V,'fro')^2;
            if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
        end
    if cnt == 1
        eps0 = nodelta; 
    end
    eps = nodelta; 
    cnt = cnt + 1; 
end