function restult =  Prox_MCP(uk, ka, lambda_pena, lambda_mcp,theta)

t     = 1/(ka*lambda_pena);
n     = length(uk);
absuk = abs(uk);
% ==============compute z========================
myfun = @(w) 0.5*(w - absuk).^2 + (lambda_mcp/t).*w - w.^2/(2*theta);
temp1 = zeros(n,1);
temp2 = theta*lambda_mcp*ones(n,1);
tfun1 = myfun(temp1);
tfun2 = myfun(temp2);
if theta == 1
    z     = temp2;
    I1    = find(tfun1 <= tfun2);
    z(I1) = temp1(I1);
else
    ttemp   = (theta/(t*theta-t)).*(t*absuk - lambda_mcp);
    temp3   = min(theta*lambda_mcp, max(0, ttemp));
    tfun3   =  myfun(temp3);
    min_fun = min(tfun1, min(tfun2, tfun3));
    
    z       = temp3;
    I1 = find(tfun1 == min_fun);
    I2 = find(tfun2 == min_fun);
    z(I1) = temp1(I1);
    z(I2) = temp1(I2);
    
end
% construct the solution
x1 = sign(uk).*z;
x2 = sign(uk).*max(theta*lambda_mcp, absuk);

mcp_fun1 = 0.5*(x1 - uk ).^2 + (1/t).* Fun_MCP(x1,lambda_mcp,theta);
mcp_fun2 = 0.5*(x2 - uk).^2  + (1/t).* Fun_MCP(x2,lambda_mcp,theta);

min_mcp_fun   = min(mcp_fun1, mcp_fun2);
Idex1         = find(mcp_fun1 == min_mcp_fun);
restult       = x2;
restult(Idex1) = x1(Idex1);
end