function z  = Fun_LSP(w,lambda_lsp,theta)
 z = lambda_lsp*log(1+abs(w)./theta);
end