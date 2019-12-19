function [tout,lh] = reduceLength(t,VaL,MaxL)
% Readme:
% If the length of input variables is less than or equal to MaxL, then the
% output variables are themselves. Otherwise, this codes helps us to adjust the
% the length of prime varibles by compute the average of some compenents.


L = length(VaL);

if L > MaxL
    I     = floor(L/MaxL);
    tout  = zeros(MaxL+1,1);
    lh    = zeros(MaxL+1,1);
    lh(1) = VaL(1);
    for j = 1 : (MaxL-1)
        tout(j+1) = mean(t(2+I*(j-1):1+I*j));
        lh(j+1)   = mean(VaL(2+I*(j-1):1+I*j));
    end
    tout(MaxL+1) = mean(t(2+I*(MaxL-1):end));
    lh(MaxL+1)   = mean(VaL(2+I*(MaxL-1):end));
else
    tout = t;
    lh   = VaL;
end
