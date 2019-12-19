function [output, lambda_pena, Num_Grad_index] = SPDCA_Table(n,X,x0,Max_grad_num,Max_prox_num,Prox)
%=========================================
% Readme:
% Optimization model
% min -sum_{i=1}^{n} \langle omega,x_i\rangle^2 + h(w) +g(w)
%where h(w) is the indicator function of {omega: sum_{j=1}^{d} omega_j = 1, omega >= 0}
% g(omega) is the sparisty function, such as LSP, Capped L1, MSP, SCAD.
%============================================
% n: the number of train set
% d: the feature number of train set
% Ka: the regularization parameter of optimization problem
% yX: y cdot X, where y is the label and X is the dataset
% x0: the initial point
%==================================================
%%% parameters which relates to data
L    = max(sum(X.^2,2));
%==========================================
% parameter of PDCA
b   = ceil(n^(1/2));
q   = b;
%==========================================
% Compute the maximum iteration number such that the total proximal and
% gradient number does not exceed the given numbers.
R   = 0;  N_grad   = 0;
while R < Max_grad_num
    if mod(N_grad,q)== 0
        R = R + n;
    else
        R = R + b;
    end
    N_grad = N_grad + 1;
end
N = min(N_grad,Max_prox_num);
%==============================================
lambda_pena = N^(-1/3);
gamma       = lambda_pena/(2*L*lambda_pena+1);
% ==============================================
% Initial point setting
TT          = [];
w           = [];
Num_Grad_index = [];
w(:,1)      = x0;
TT(1)       = 0;
Num_Grad_index(1) = 0;
wk_old_old  = x0;
wk_old      = x0;
%===============================================
num_grad = 0;
num_prox = 0;
ErrRun  = 1;
k       = 1;
tic
while ErrRun == 1
    if mod(k-1,q) == 0
        num_grad        = num_grad +  n;
        Num_Grad_index  = [Num_Grad_index num_grad];
        arg             =  X*wk_old;
        Xarg            = (X.*arg);
        grad_esti       = -(sum(Xarg,1)/n)';
    else
        num_grad         = num_grad +  b;
        Num_Grad_index   = [Num_Grad_index num_grad];
        
        I          =  randi(n,[b,1]);
        XI         = X(I,:);
        argI1      =  XI*wk_old;
        argI2      =  XI*wk_old_old;
        tempI1     =  (XI.*argI1);
        tempI2     =  (XI.*argI2)  ;
        
        Diff_grad  = (-sum(tempI1,1)/b)' +  (sum(tempI2,1)/b)';
        grad_esti  = grad_esti +  Diff_grad;
    end
    subg     = Prox(wk_old,lambda_pena);
    prox_vec = wk_old - gamma*(grad_esti+ (wk_old - subg)./lambda_pena);
    wk_new   = BaProj(prox_vec);
    
    num_prox  = num_prox + 1;
    if (num_prox > Max_prox_num) || (num_grad > Max_grad_num)
        ErrRun = 0;
    end
    fprintf('--SPDCA--num_grad/Max_grad_num=%1.4f---num_prox/Max_prox_num=%1.4f---\n',num_grad/Max_grad_num, num_prox/Max_prox_num );
    
    wk_old_old = wk_old;
    wk_old     = wk_new;
    
    k   = k + 1;
    w   = [w wk_new];
    TT  = [TT toc];
end
output = [w; TT];
end