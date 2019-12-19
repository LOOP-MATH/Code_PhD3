function [output,lambda_pena] = VRSPA(n,X,x0,IsGradCom,Ngpc,Prox)
%=========================================
% Readme:
%model
% min - sum_{i=1}^{n} \langle omega,x_i\rangle^2 + h(w) +g(w)
%where h(w) is the indicator function of {omega: sum_{j=1}^{d} omega_j = 1, omega >= 0}
% g(omega) is the sparisty function which satisfy g_1 - g_2.
% Reference information:
% For more details, see
% Metel and Takeda. Stochastic gradient methods for non-smooth nonconvex
% reguralized optimization. ICML 2019.
%============================================
% n: the number of train set
% d: the feature number of train set
% Ka: the regularization parameter of optimization problem
% yX: y cdot X, where y is the label and X is the dataset
% x0: the initial point
% IsGradCom: If IsGradCom==1, for gradient complexity, IsGradCom =0, for
% proximal complexity
% Ngpc: the index which controls the maximum number of computation of gradient or proximal operator.
% L1: for compute the Lipschitz constant of the smooth function
%=========================================
m      = ceil(n^(1/3));
b      = m^2;
if IsGradCom == 1
    N = ceil(Ngpc/(n/m+b));
else
    N = Ngpc;
end
%=========================================
S      = ceil(N/m);
lambda_pena = (S*m)^(-1/3);

L      = max(sum(X.^2,2));
LEr    = L+1/lambda_pena;
beta   = 1/(6*LEr);
%=========================================
% setting the initial point
TT     = [];
w      = [];

w(:,1)  = x0;
TT(1)   = 0;
num     = 1;
ErrRun  = 1;
k       = 1;
%=========================================
tic
while ErrRun == 1
    w1  = w(:,1+m*(k-1));
    arg =  X*w1;
    Gd  =  (X.*arg);%Gradient decomposed
    G   = (-sum(Gd,1)/n)';
    if IsGradCom==1
        num  = num +  n;
        if num > 0.995*Ngpc
            ErrRun = 0;
        end
        fprintf('-----VRSPA------num/Ngpc=%1.4f-------\n',num/Ngpc);
    end
    t = 1;
    while t <= m && (ErrRun==1)
        I         = randi(n,[b,1]);
        wkt       = w(:,t+m*(k-1));
        
        arg       = X(I,:)*wkt;
        GDI       = (X(I,:).*arg);
        grad_esti = -(sum(GDI,1)/b)' + (sum(Gd(I,:),1)/b)' + G;
        
        xi        = Prox(wkt,lambda_pena);
        vrg       = (wkt - xi)/lambda_pena + grad_esti;
        if IsGradCom == 1
            num  = num +  b;
            if num > 0.995*Ngpc
                ErrRun = 0;
            end
        else
            num  = num + 1;
            if num > Ngpc
                ErrRun = 0;
            end
        end
        fprintf('-----VRSPA------num/Ngpc=%1.4f-------\n',num/Ngpc);
        
        wk_new = BaProj(wkt-beta*vrg);
        w   = [w wk_new];
        TT  = [TT toc];
        
        t = t + 1;
    end
    k = k + 1;
end
%=========================================
output = [w; TT];
end