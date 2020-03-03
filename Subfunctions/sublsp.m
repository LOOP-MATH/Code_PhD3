function z = sublsp(x, ka, lambda, theta)
temp1 = subL1(x,1,1);
temp2 = (1./(1+ abs(x)/theta)).*(temp1/theta);
z     = ka*lambda*(temp1 - temp2);
end