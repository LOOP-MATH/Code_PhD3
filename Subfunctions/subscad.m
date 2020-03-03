function z = subscad(x, ka, lambda, theta)

tem_abs = abs(x);
id1    = find(tem_abs <= lambda);
id3    = find(tem_abs > theta*lambda);

temp3 = subL1(x,1,lambda);
y     = (2*x - 2*temp3)./((theta-1)*2);
 
y(id1,:) = zeros(length(id1),1);
y(id3,:) = temp3(id3,:);

z = ka*y;
end