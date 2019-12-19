function x = Prox_CappedL1(uk,ka, lambda_pena,lambda_capp,theta)
% The formulation was given by  Lu et al.
% min_x .5||x-y||^2 + lambda/t sum_i min(\|w_i\|, theta)
% lambda is the parameter of CappedL1.
% t is the penalty parameter.

t     = 1/(ka*lambda_pena);
absuk = abs(uk);

x1    = sign(uk).*max(theta, absuk);
x2    = sign(uk).*min(theta, max( absuk - lambda_capp/t, 0));

h1   = 1/2*(x1-uk).^2 + (1/t).*Fun_CappL1(x1,lambda_capp,theta);
h2   = 1/2*(x2-uk).^2 + (1/t).*Fun_CappL1(x2,lambda_capp,theta);

x    = x1;
I    = find(h1 > h2);
x(I) = x2(I);
end