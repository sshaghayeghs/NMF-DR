function [W,H,e,t] = PGLIN(V,Winit,Hinit,tol,timelimit,maxiter)

% NMF by alternative non-negative least squares using projected gradients
% Author: Chih-Jen Lin, National Taiwan University
% Source: Lin, Projected Gradient Methods for Nonnegative Matrix Factorization,
% Neural Computation, 19, p. 2756-2779, 2007, MIT press. 
% W,H: output solution
% Winit,Hinit: initial solution
% tol: tolerance for a relative stopping condition
% timelimit, maxiter: limit of time and iterations
initt = cputime;  nM = norm(V,'fro')^2; 
W = Winit; H = Hinit; t = []; e = [];
% We also scale initial iterates, otherwise does not work well
HHT = H*H'; WTW = W'*W; VHT = V*H'; 
gradW = W*(HHT) - VHT; gradH = (WTW)*H - W'*V;
initgrad = norm([gradW; gradH'],'fro');
%fprintf('Init gradient norm %f\n', initgrad);
tolW = max(0.001,tol)*initgrad; tolH = tolW;
for iter=1:maxiter,
    % stopping condition
    projnorm = norm([gradW(gradW<0 | W>0); gradH(gradH<0 | H>0)]);
    if projnorm < tol*initgrad | cputime-initt > timelimit,
        break;
    end
    [W,gradW,iterW,HHt,VtH] = nlssubprob(V',H',W',tolW,1000,cputime-initt);
    cnT = cputime; 
    e = [e sqrt( (nM-2*sum(sum(W.*(VtH)))+ sum(sum(HHt.*(W*W')))) )];
    initt = initt+(cputime-cnT);
    t = [t cputime-initt]; 
    W = W'; gradW = gradW'; 
    if iterW==1,
        tolW = 0.1 * tolW;
    end
    [H,gradH,iterH,WtW,WtV] = nlssubprob(V,W,H,tolH,1000,cputime-initt);
    cnT = cputime; 
    e = [e sqrt( (nM-2*sum(sum(H.*(WtV)))+ sum(sum(WtW.*(H*H')))) )];
    initt = initt+(cputime-cnT);
    t = [t cputime-initt]; 
    if iterH==1,
        tolH = 0.1 * tolH;
    end
    if rem(iter,10)==0, fprintf('.'); end
end
%fprintf('\nIter = %d Final proj-grad norm %f\n', iter, projnorm);

function [H,grad,iter,WtW,WtV] = nlssubprob(V,W,Hinit,tol,maxiter,tim)
% H, grad: output solution and gradient
% iter: #iterations used
% V, W: constant matrices
% Hinit: initial solution
% tol: stopping tolerance
% maxiter: limit of iterations
initt = cputime; 
H = Hinit; WtV = W'*V; WtW = W'*W; grad = WtW*H - WtV; 
alpha = 1; beta = 0.1; iter = 1;
while (iter == 1) || (iter <= maxiter && cputime-initt < tim),
    grad = WtW*H - WtV;
    projgrad = norm(grad(grad < 0 | H >0));
    if projgrad < tol 
        break
    end
    % search step size
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
        if cputime-initt > tim,
            break
        end
    end
    iter = iter + 1;
end