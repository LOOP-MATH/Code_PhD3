function z  = Fun_MCP(w,lambda_mcp,theta)
absw   = abs(w);
Id1    = find(absw <= theta*lambda_mcp);
z      =  (theta*lambda_mcp^2/2)*ones(length(w),1);
tmp1   = lambda_mcp.*abs(w) - w.^2/(2*theta);
z(Id1) = tmp1(Id1);
end