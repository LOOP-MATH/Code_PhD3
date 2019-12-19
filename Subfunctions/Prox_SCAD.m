function z =  Prox_SCAD(uk, ka, lambda_pena, lambda_scad,theta)

t     = 1/(ka*lambda_pena);
absuk = abs(uk);

x1    = sign(uk).*min(lambda_scad, max(0, absuk-lambda_scad/t ));
x3    = sign(uk).*max(theta*lambda_scad, absuk);
temp1 = (t*(theta-1)*absuk - theta*lambda_scad)/(t*(theta-2));
x2    = sign(uk).*min(theta*lambda_scad, temp1);

h1 = 0.5*norm(x1-uk,2).^2 + 1/t.* Fun_SCAD(x1,lambda_scad,theta);
h2 = 0.5*norm(x2-uk,2).^2 + 1/t.* Fun_SCAD(x2,lambda_scad,theta);
h3 = 0.5*norm(x3-uk,2).^2 + 1/t.* Fun_SCAD(x3,lambda_scad,theta);

z     = zeros(length(uk),1);
Temp2 = min(h1,min(h2,h3));
I1    = find(h1 == Temp2);
z(I1) = x1(I1);
I2    = find(h2==Temp2);
z(I2) = x2(I2);
I3    = find(h3==Temp2);
z(I3) = x1(I3);
end