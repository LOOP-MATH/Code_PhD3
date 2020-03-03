function output = SSPDCA_Figure(n,X,x0,subgradient_fun_f1,...
    subgradient_fun_f2,Max_smoothgrad_num,Max_nonsmoothsub_num)
%===========================================
% Readme:
% SSDC + SPG
% min  (f_1(w) - f_2(w)) +  fsmooth(w), s.t. w\in K
% where f_1 f_2 are smooth convex functions,
% K is a nonempty closed convex set.
% ======= Decomposition 1======
% Relaxed resulting subproblem I:
% rho* f_1(w) - \langle\xi_k, w\rangle + fsmooth(w) + gamma/2 \|w-wk\|^2
% For more details, see
% Xu et al. Stochastic optimization for DC Functions and non-smooth
% non-convex regularizers with non-asymptotic convergence.
%===========================================
L     = max(sum(X.^2,2));
alpha = 0.99/((1+2*sqrt(2))*L);
%===========================================
sum_grad_smooth    = 0;
sum_subg_nonsmooth = 0;
%==========================================
% setting the initial point
TT          = [];
w           = [];
w(:,1)      = x0;
TT(1)       = 0;
%==========================================
ErrRun     = 1;
wk_old     =  x0;
wk_old_old = x0;
q = floor(sqrt(n));
b = q;


Index_sum = 1;
k       = 1;
tic
while ErrRun == 1
    subgradientf2_xi = subgradient_fun_f2(wk_old);
    if mod(k-1,q) == 0
        sum_grad_smooth = sum_grad_smooth +  n;
        arg             =  X*wk_old;
        Xarg            = (X.*arg);
        grad_esti       = -(sum(Xarg,1)/n)';
    else
        sum_grad_smooth = sum_grad_smooth +  b;
        
        I          =  randi(n,[b,1]);
        XI         = X(I,:);
        argI1      =  XI*wk_old;
        argI2      =  XI*wk_old_old;
        tempI1     =  (XI.*argI1);
        tempI2     =  (XI.*argI2)  ;
        
        Diff_grad_smooth  = (-sum(tempI1,1)/b)' +  (sum(tempI2,1)/b)';
        grad_esti         = grad_esti +  Diff_grad_smooth;
    end
    pk = wk_old - alpha*(grad_esti- subgradientf2_xi);
    if sum_grad_smooth > Max_smoothgrad_num
        ErrRun = 0 ;
    end
    %%
    % compute the subproblem.
    t         = 1;
    Tk        = k +1;  %2
    while t <= Tk && ErrRun == 1
        eta   = 1/(t+1);
        if t== 1
            wkt   = wk_old;
        end
        prox_vec = subgradient_fun_f1(wkt)+ 1/alpha*(wkt-pk);
        wkt_new  = BaProj(wkt - eta*prox_vec);
        t        = t+1;
        wkt      = wkt_new;
        
        w   = [w wkt_new];
        TT  = [TT toc];
        sum_subg_nonsmooth = sum_subg_nonsmooth +1;
        if sum_grad_smooth > Max_smoothgrad_num ||  sum_subg_nonsmooth  > Max_nonsmoothsub_num
            ErrRun = 0 ;
        end
        fprintf('--SSPDCA--num_Smooth_grad/Max_grad_num=%1.4f---num_nonsmooth_gra/Max_nonsmooth_num=%1.4f---\n',...
            sum_grad_smooth/Max_smoothgrad_num, sum_subg_nonsmooth/Max_nonsmoothsub_num );
    end
    if ErrRun == 1
        index     = Index_sum:Index_sum + Tk-1;
        index1    = 1:Tk;
        wk_new    = sum(w(:,index).* (index1),2)/( sum(index1));
        Index_sum = Index_sum + Tk;
        w(:, Index_sum) = wk_new;
        
        wk_old_old  = wk_old;
        wk_old      = wk_new;
        
        k = k + 1;
    end
end
output = [w;TT];
end