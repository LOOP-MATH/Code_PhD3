function [output, Num_Grad_index,Num_subGrad_index] = SSPDCA_Table(n,X,x0,...
    Max_grad_num,Max_subgradient_num,...
    subgradient_fun_f2,kalambda)
%===========================================
% Readme:
%===========================================
L             = max(sum(X.^2,2));
alpha         = 1/(2*L);
b             = ceil(sqrt(n));
q             = b;
kalambdaalpha = alpha*kalambda;
ex            = ones(length(x0),1);
%==========================================
% setting the initial point
TT          = [];
w           = [];
Num_Grad_index = [];
Num_Grad_index(1) = 0;
Num_subGrad_index = [];
Num_subGrad_index(1) = 0;

w(:,1)      = x0;
TT(1)       = 0;


wk_old      = x0;
wk_old_old  = x0;
%==========================================
num_grad   = 0;
num_subg   = 0;
ErrRun     = 1;
k          = 1;
tic
while ErrRun == 1
    if mod(k-1,q) == 0
        num_grad = num_grad +  n;
        Num_Grad_index  = [Num_Grad_index num_grad];
        arg             =  X*wk_old;
        Xarg            = (X.*arg);
        grad_estimation = -(sum(Xarg,1)/n)';
    else
        num_grad        = num_grad +  b;
        Num_Grad_index  = [Num_Grad_index num_grad];
        
        I          =  randi(n,[b,1]);
        XI         = X(I,:);
        argI1      =  XI*wk_old;
        argI2      =  XI*wk_old_old;
        tempI1     =  (XI.*argI1);
        tempI2     =  (XI.*argI2)  ;
        
        Diff_grad_smooth  = (-sum(tempI1,1)/b)' +  (sum(tempI2,1)/b)';
        grad_estimation   = grad_estimation +  Diff_grad_smooth;
    end
    
    subgradientf2_xi = subgradient_fun_f2(wk_old);
    
    Pk      = wk_old - alpha*(grad_estimation- subgradientf2_xi);
    Proj    = Pk  - kalambdaalpha*ex;
    wk_new  = BaProj(Proj);
    
    num_subg = num_subg +1;
    Num_subGrad_index  = [Num_subGrad_index num_subg];
    if (num_subg > Max_subgradient_num) || (num_grad > Max_grad_num)
        ErrRun = 0;
    end
    fprintf('--SSPDCA--num_Smooth_grad/Max_grad_num=%1.4f---num_nonsmooth_gra/Max_nonsmooth_num=%1.4f---\n',...
        num_grad/Max_grad_num, num_subg/Max_subgradient_num );
    
    wk_old_old  = wk_old;
    wk_old      = wk_new;
    
    
    k   = k + 1;
    w   = [w wk_new];
    TT  = [TT toc];
    
end
output = [w; TT];
end