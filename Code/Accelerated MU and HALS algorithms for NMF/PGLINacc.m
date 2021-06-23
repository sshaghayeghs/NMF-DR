function [W,H,e,t] = PGLINacc(V,Winit,Hinit,alpha,delta,maxiter,timelimit)

% Modified version of the 
%
% NMF by alternative non-negative least squares using projected gradients
% Author: Chih-Jen Lin, National Taiwan University
% Source: Lin, Projected Gradient Methods for Nonnegative Matrix Factorization,
% Neural Computation, 19, p. 2756-2779, 2007, MIT press. 
% W,H: output solution
% Winit,Hinit: initial solution
% alpha: parameter for acceleration
% delta: parameter for stopping inner iterations
% timelimit, maxiter: limit of time and iterations
%
% using a fixed number of inner iterations (depending of parameter alpha)

initt = cputime;
nM = norm(V,'fro')^2; [m,n] = size(V); [m,r] = size(Winit);

if nargin <= 3, alpha = 0.5; end
if nargin <= 4, delta = 0; end
if nargin <= 5, maxiter = 100; end
if nargin <= 6, timelimit = 60; end

if issparse(V), K = sum(V(:) ~= 0); else K = m*n; end 
rhoW = 1+floor((K+n*r)/(m*(r+1)));
rhoH = 1+floor((K+m*r)/(n*(r+1)));
W = Winit; H = Hinit;  t = []; e = []; iter = 0;

while iter <= maxiter && cputime-initt < timelimit,
    [W,gradW,iterW,HHt,VtH] = nlssubprob(V',H',W',1+alpha*rhoW,delta); 
    cnT = cputime; 
    e = [e sqrt( (nM-2*sum(sum(W.*(VtH)))+ sum(sum(HHt.*(W*W')))) )];
    initt = initt+(cputime-cnT);
    t = [t cputime-initt]; 
    W = W'; 
    
    [H,gradH,iterH,WtW,WtV] = nlssubprob(V,W,H,1+alpha*rhoH,delta);
    cnT = cputime; 
    e = [e sqrt( (nM-2*sum(sum(H.*(WtV)))+ sum(sum(WtW.*(H*H')))) )];
    initt = initt+(cputime-cnT);
    t = [t cputime-initt]; 
    iter = iter + 1;
end

function [H,grad,iter,WtW,WtV] = nlssubprob(V,W,Hinit,maxiter,delta)
% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations
H = Hinit; WtV = W'*V; WtW = W'*W; 
% At first pass, scale solution (p. 72 of the thesis)
alpha = 1; beta = 0.1; iter = 1; eps = 1; eps0 = 1;
while iter <= maxiter && eps >= delta*eps0,
    grad = WtW*H - WtV;
    projgrad = norm(grad(grad < 0 | H >0));
    % search step size
    H0 = H; 
    for inner_iter=1:20,
        Hn = max(H - alpha*grad, 0); d = Hn-H;
        gradd=sum(sum(grad.*d)); dQd = sum(sum((WtW*d).*d));
        suff_decr = 0.99*gradd + 0.5*dQd < 0;
        if inner_iter==1,
            decr_alpha = ~suff_decr; Hp = H;
        end
        if decr_alpha,
            if suff_decr,
                H = Hn; break;
            else
                alpha = alpha * beta;
            end
        else
            if ~suff_decr | Hp == Hn,
                H = Hp; break;
            else
                alpha = alpha/beta; Hp = Hn;
            end
        end
    end
    if iter == 1
        eps0 = norm(H-H0,'fro');
    end
    eps = norm(H-H0,'fro');
    iter = iter+1; 
end