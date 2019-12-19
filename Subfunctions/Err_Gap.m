function zGap = Err_Gap(w,X,lambda_pena,Prox,Mygfun,Name)
n  = size(X,1);
it = size(w,2);
zGap  = zeros(it,1);
for k = 1:it
    wk   = w(:,k);
 
    xi       =  Prox(wk, lambda_pena);
    arg      =  X*xi;
    Gd       =  (X.*arg);%Gradient decomposed
    grad_xi  = -(sum(Gd,1)/n)';
    proj1    = xi - (grad_xi+ (wk-xi)/lambda_pena);
    tildewk  = BaProj(proj1);
    zGap(k)     = norm(wk - tildewk);
    
    disp(['---',Name,'---iter = ', num2str(k/it),'---ErrVal = ', num2str(zGap(k)),'---']); 
end
end