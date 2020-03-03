function zFun = Err_Fun_Gap_SS(w,X,Mygfun,Name)
n  = size(X,1);
it = size(w,2);
zFun  = zeros(it,1);
 
for k = 1:it
    wk   = w(:,k);
    arg  =  (X*wk);
    f    = -sum(arg.*arg)/(2*n);
    g    =   Mygfun(wk);
    zFun(k) = (f+g);
    disp(['---',Name,'---iter = ', num2str(k/it),'---FunVal = ', num2str(zFun(k)),'---']);
end
end