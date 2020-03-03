function z = subl1l2(x, ka, lambda_l1l2, theta)
temp = norm(x);
if temp < 1e-6
    z = zeros(length(x),1);
else
    z =  (ka*lambda_l1l2*theta/temp)*x;
end
end