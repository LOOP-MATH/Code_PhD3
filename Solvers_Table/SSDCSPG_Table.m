function [output, Num_Grad_index, Num_subGrad_index]  = SSDCSPG_Table(n,X,x0,b1,subgradient_fun_f1,...
    subgradient_fun_f2,Max_grad_num,Max_subgradient_num)
%===========================================
% Readme:
% SSDC + SPG
% min rho*(f_1(w) - f_2(w)) +  fsmooth(w), s.t. w\in K
% where f_1 f_2 are smooth convex functions,
% K is a nonempty closed convex set.
% ======= Decomposition 1======
% Relaxed resulting subproblem I:
% rho* f_1(w) - \langle\xi_k, w\rangle + fsmooth(w) + gamma/2 \|w-wk\|^2
% For more details, see
% Xu et al. Stochastic optimization for DC Functions and non-smooth
% non-convex regularizers with non-asymptotic convergence.
%===========================================
L      = max(sum(X.^2,2));
gamma  = 3*L;
%==========================================
% setting the initial point
TT     = [];
w      = [];
Num_Grad_index = [];
Num_Grad_index(1) = 0;

Num_subGrad_index = [];
Num_subGrad_index(1) = 0;

w(:,1) = x0;
TT(1)  = 0;
%==========================================
num_grad = 0;
num_subg = 0;
ErrRun = 1;
k       = 1;
Index_sum = 1;
tic
while ErrRun == 1
    
    w1        =  w(:, Index_sum);
    subf2_xi  =  subgradient_fun_f2(w1);
    
    t         = 1;
    Tk        = ceil(k/3)+1;
    while t <= Tk && ErrRun == 1
        etat       = 4/(gamma*t);
        gamma_etat = gamma + 1/etat;
        wkt        = w(:,Index_sum+t-1);
        %----------------------------
        % subgradient method for subproblem
        I        = randi(n,[b1,1]);
        XI       = X(I,:);
        argI1    =  XI*wkt;
        tempI1   =  (XI.*argI1);
        prox_vec = (-sum(tempI1,1)/b1)'+ subgradient_fun_f1(wkt)- subf2_xi - gamma*w1 - 1/etat*wkt;
        wkt_new  = BaProj(-prox_vec/gamma_etat);
        num_subg = num_subg +1;
        Num_subGrad_index  = [Num_subGrad_index num_subg];
        
        w        = [w wkt_new];
        TT       = [TT toc];
        t        = t+1;
        
        num_grad = num_grad + b1;
        Num_Grad_index  = [Num_Grad_index num_grad];
        num_subg = num_subg + 1;
        if num_grad > Max_grad_num ||  num_subg  > Max_subgradient_num
            ErrRun = 0 ;
        end
        fprintf('--SSDCSPG--num_Smooth_grad/Max_grad_num=%1.4f---num_nonsmooth_gra/Max_nonsmooth_num=%1.4f---\n',...
            num_grad/Max_grad_num, num_subg/Max_subgradient_num );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ErrRun == 1
        index     = Index_sum:Index_sum + Tk-1;
        idex1     = 1:Tk;
        wk_new    = sum(w(:,index).* (idex1),2)/sum(idex1);
        w(:,end)  = wk_new;
        Index_sum = Index_sum + Tk;
        k = k+1;
    end
    
    
end
output = [w;TT];
end