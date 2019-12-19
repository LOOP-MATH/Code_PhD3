function x = Prox_L1L2(uk,ka,lambda_pena,lambda_l1l2,theta)
% download from  Ming Yan's paper.
% min_x .5||x-y||^2 + lambda(|x|_1- theta |x|_2)

x = zeros(length(uk));
t = ka*lambda_pena*lambda_l1l2;

if max(abs(uk)) > 0
    if max(abs(uk)) > t
        x   = max(abs(uk) - t, 0).*sign(uk);
        x   = x * (norm(x) + theta * t)/norm(x);
    else
        if max(abs(uk))>=(1-theta)*t
            [~, i]  = max(abs(uk));
            x(i(1)) = (abs(uk(i(1))) + (theta - 1) * t) * sign(uk(i(1)));
        end
    end
end
end