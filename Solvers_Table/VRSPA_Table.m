function [output, lambda_pena,Num_Grad_index]  = VRSPA_Table(n,X,x0,Max_grad_num,Max_prox_num,Prox)
%=========================================
% Readme:
% Optimization model
% min -sum_{i=1}^{n} \langle omega,x_i\rangle^2 + h(w) +g(w)
%where h(w) is the indicator function of {omega: sum_{j=1}^{d} omega_j = 1, omega >= 0}
% g(omega) is the sparisty function, such as LSP, Capped L1, MSP, SCAD.
% Reference information:
% For more details, see
% Metel and Takeda. Stochastic gradient methods for non-smooth nonconvex
% reguralized optimization. 2019.
%============================================
% n: the number of train set
% d: the feature number of train set
% Ka: the regularization parameter of optimization problem
% yX: y cdot X, where y is the label and X is the dataset
% x0: the initial point
%=========================================
m      = ceil(n^(1/3));
b      = m^2;
% Compute the maximum iteration number such that the total proximal and
% gradient number does not exceed the given numbers.
N_grad = ceil(Max_grad_num/(n/m+b));
N      = min(N_grad,Max_prox_num);
%=========================================
S           = ceil(N/m);
lambda_pena = (S*m)^(-1/3);

L      = max(sum(X.^2,2));
LEr    = L+1/lambda_pena;
beta   = 1/(6*LEr);
%=========================================
% setting the initial point
TT     = [];
w      = [];
Num_Grad_index = [];
w(:,1)   = x0;
TT(1)    = 0;
Num_Grad_index(1) = 0;
num_grad = 0;
num_prox = 0;
ErrRun   = 1;
k        = 1;
%=========================================
tic
while ErrRun == 1
    w1  = w(:,1+m*(k-1));
    arg =  X*w1;
    Gd  =  (X.*arg);
    G   = (-sum(Gd,1)/n)';
    num_grad  = num_grad +  n;
    if num_grad > Max_grad_num
        ErrRun = 0;
    end
    t = 1;
    while t <= m && (ErrRun==1)
        I         = randi(n,[b,1]);
        num_grad  = num_grad +  b;
        Num_Grad_index   = [Num_Grad_index num_grad];
        num_prox  = num_prox + 1;
        if ( num_grad > Max_grad_num) || ( num_prox > Max_prox_num)
            ErrRun = 0;
        end
        fprintf('--VRSPA--num_grad/Max_grad_num=%1.4f---num_prox/Max_prox_num=%1.4f---\n',num_grad/Max_grad_num, num_prox/Max_prox_num );
        
        
        wkt       = w(:,t+m*(k-1));
        xi        = Prox(wkt,lambda_pena);
        
        arg       = X(I,:)*wkt;
        Grad_cur  = (X(I,:).*arg);
        grad_esti = -(sum(Grad_cur,1)/b)' + (sum(Gd(I,:),1)/b)' + G;
        vrg       =  grad_esti + (wkt - xi)/lambda_pena;
        
        wk_new = BaProj(wkt-beta*vrg);
        w   = [w wk_new];
        TT  = [TT toc];
        
        t = t + 1;
    end
    k = k + 1;
end
output = [w; TT];
end