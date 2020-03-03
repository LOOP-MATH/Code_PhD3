function [w,sum_subg_nonsmooth] = StochasticSubgrad(y0, Tk,Pk,mu,kalambda,alpha,...
    sum_subg_nonsmooth,Max_nonsmooth_subgradient_num)

y = [];
y_old = y0;
y = [y y_old];
t = 0;
while t <= Tk && sum_subg_nonsmooth <= Max_nonsmooth_subgradient_num
    idex= find( abs(y_old)<= 1e-6 );
    temp = kalambda*sign(y_old);
    temp(idex,:) = zeros(length(idex),1);
    
    xi1  = temp + 1/alpha *(y_old - Pk);
    gammat = 2/mu/(t+1);
    y_new  = BaProj(y_old- gammat*xi1);
    y_old  = y_new;
    y = [y y_old];
    t = t +1;
    sum_subg_nonsmooth = sum_subg_nonsmooth +1;
end

Y_index = y(:,1:end-1);
[~,s] =size(Y_index);
Index1 = 1: s;

z =  Y_index.*Index1;

w = 2*sum(z,2)./(s*(s+1));

end