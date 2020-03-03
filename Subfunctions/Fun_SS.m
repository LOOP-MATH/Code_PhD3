function zFun = Fun_SS(n,w,X,Mygfun)
    arg  =  (X*w);
    f    = -sum(arg.*arg)/(2*n);
    g    =   Mygfun(w);
    zFun  = (f+g);
end