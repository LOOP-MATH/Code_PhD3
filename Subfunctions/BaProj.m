function  z =  BaProj(w)
temp = max(w,0);
c    = max(norm(temp),1);
z    = (1/c).* temp;
end