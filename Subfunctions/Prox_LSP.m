function restult =  Prox_LSP(uk, ka, lambda_pena, lambda_lsp,theta)
t     = 1/(ka*lambda_pena);
absuk = abs(uk);
n     = length(uk);

proxfunval = @(w) 0.5*(w- absuk).^2 + (lambda_lsp/t)*log(1+w./theta);

% compute z
x1     = zeros(n,1);
temp0  =  t^2.*(absuk - theta).^2 - 4*t*(lambda_lsp - (t*theta).*absuk);
Index1 = find(temp0 >= 0);
tukminustheta =  t.*(absuk-theta);
x2         = x1;
x2(Index1) = max((tukminustheta(Index1) + sqrt(temp0(Index1)))./(2*t),0);
x3 = x1;
x3(Index1) = max(tukminustheta(Index1) - sqrt(temp0(Index1))./(2*t), 0);

fun1 = proxfunval(x1);
fun2 = proxfunval(x2);
fun3 = proxfunval(x3);
min_fun = min(fun1, min(fun2,fun3));
Idex2   = find(min_fun == fun2);
Idex3   = find(min_fun == fun3);

restult = x1;
restult(Idex2)= x2(Idex2);
restult(Idex3)= x2(Idex3);
end