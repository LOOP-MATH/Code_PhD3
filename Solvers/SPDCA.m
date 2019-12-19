function [output, lambda_pena] = SPDCA(n,X,x0,IsGradCom,Ngpc,Prox)
% n: the number of train set
% d: the feature number of train set
% Ka: the regularization parameter of optimization problem
% X:  train data set
% x0: the initial point
% IsGradCom: If IsGradCom==1, for gradient complexity, IsGradCom =0, for
% proximal complexity
% Ngpc: the index which controls the maximum number of computation of gradient or proximal operator.
% theta
%==================================================
%%% parameters which relates to data
L    = max(sum(X.^2,2));
%==========================================
% parameter of PDCA
b   = ceil(n^(1/2));
q   = b;
%==========================================
% proximal complexity and gradient complexity
if IsGradCom == 1
    R   = 0;   %actual number of gradient calls used, which is close to and >= N
    i   = 0;
    while R < Ngpc
        if mod(i,q)== 0
            R = R + n;
        else
            R = R + b;
        end
        i = i+1;
    end
    N = i;
else
    N =  Ngpc;
end
%==============================================
lambda_pena = N^(-1/3);
gamma  = lambda_pena/(2*L*lambda_pena+1);
% ==============================================
% Initial point setting
TT         = [];
w          = [];
w(:,1)     = x0;
TT(1)      = 0;
wk_old_old = x0;
wk_old     = x0;
%===============================================

num     = 1;
ErrRun  = 1;
k       = 1;
tic
while ErrRun == 1
    if mod(k-1,q) == 0
        arg      =  X*wk_old;
        Xarg    = (X.*arg);
        vk       = -(sum(Xarg,1)/n)';
        if IsGradCom==1
            num  = num +  n;
            if num > 0.995*Ngpc
                ErrRun = 0;
            end
            fprintf('-----SPDCA------num/Ngpc=%1.4f-------\n',num/Ngpc);
        end
    else
        if IsGradCom==1
            num  = num +  b;
            if num > 0.995*Ngpc
                ErrRun = 0;
            end
            fprintf('-----SPDCA------num/Ngpc=%1.4f-------\n',num/Ngpc);
        end 
        I          =  randi(n,[b,1]);
        XI         = X(I,:);
        
        argI1      =  XI*wk_old;
        argI2      =  XI*wk_old_old;
        
        tempI1     =  (XI.*argI1);
        tempI2     =  (XI.*argI2)  ;

        Diff_grad  = (-sum(tempI1,1)/b)' +  (sum(tempI2,1)/b)';
        vk         = vk +  Diff_grad;
    end
        
   
    subg     = Prox(wk_old,lambda_pena);
    
    prox_vec = wk_old - gamma*(vk+ (wk_old - subg)./lambda_pena);
    wk_new   = BaProj(prox_vec);
    if IsGradCom==0
        num  = num +  1;
        if num >  Ngpc
            ErrRun = 0;
        end
        fprintf('-----SPDCA------num/Ngpc=%1.4f-------\n',num/Ngpc);
    end
    
    
    wk_old_old = wk_old;
    wk_old     = wk_new;
    
    k   = k +1;
    w   = [w wk_new];
    TT  = [TT toc];
end
output = [w; TT];
end