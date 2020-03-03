function z  = submcp(x, ka, lambda, theta)
temp_abs = abs(x);
id2      = find(temp_abs > lambda*theta);

temp2   = subL1(x,1,lambda);
 y       = x./theta;
y(id2,:) = temp2(id2,:);
z        = ka*y;
end