function z  = Fun_SCAD(w,lambda_scad,theta)

if theta <= 2
    fprintf('theta should be great than 2.');
end

absw = abs(w);
Id1  = find(absw <= lambda_scad);
Id3  = find(absw > theta*lambda_scad);

z = (-w.^2+ 2*theta*lambda_scad.*absw - lambda_scad^2)/(2*theta-2);
z(Id1) = lambda_scad*absw(Id1);
z(Id3) = (theta+1)*lambda_scad^2/2;
end