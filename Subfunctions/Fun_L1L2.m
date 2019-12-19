function z = Fun_L1L2(wk,lambda_l1l2,theta)
z = lambda_l1l2*(norm(wk,1)-theta*norm(wk,2));
end