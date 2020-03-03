function [zFun,zGap,zTime_stop,zTime_iter,zNum_grad, zNum_prox] =  ...
    Err_Fun_Gap_Table_SS(w,TT,Num_Grad,Num_subGrad,X,lambda_pena,Prox,Mygfun,Name,eps)
%==========================================================================
n        = size(X,1);
num_iter = size(w,2);
Gap   = 1;
k     = 1;
zTime_stop = 0;
while k <= num_iter && Gap > eps
    wk      = w(:,k);
    tic
    xi       =  Prox(wk,lambda_pena);
    arg      =  X*xi;
    Gd       =  (X.*arg);    %Gradient decomposed
    grad_xi  = -(sum(Gd,1)/n)';
    proj1    = xi - (grad_xi+ (wk-xi)/lambda_pena);
    tildewk  = BaProj(proj1);
    Gap      = norm(wk - tildewk);
    zTime_stop    = zTime_stop + toc;
    disp(['---',Name,'---iter = ', num2str(k/num_iter),'---ErrVal = ', num2str(Gap),'----']);
    k        = k +1;
end
zTime_iter = TT(k-1);
arg  =  (X*wk);
f    = -sum(arg.*arg)/(2*n);
g    =   Mygfun(wk); %min(abs(wk),theta);
zFun = (f+g);
zGap = Gap;
zNum_grad = Num_Grad(k-1);
zNum_prox = Num_subGrad(k-1);
end