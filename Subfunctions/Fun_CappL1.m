function z = Fun_CappL1(wk,lambda_capp,theta)
z     = lambda_capp*min(abs(wk),theta);
end