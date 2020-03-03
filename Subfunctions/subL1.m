function z =  subL1(x,ka,lambda)
%id  = find(abs(x)<= 1e-6);
z = ka*lambda.*sign(x);
%z(id,:) = zeros(length(id),1);
end