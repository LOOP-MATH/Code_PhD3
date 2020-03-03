function z = subCapp(x, ka, lambda_capp, theta)

diff    = abs(x) - theta;
z       = zeros(length(x),1);
id      = find(diff > 0);
temp1   = sign(x);
z(id,:) = ka*lambda_capp.*temp1(id,:);

end