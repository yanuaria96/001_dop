function [R theta phi] = Cart2Sph(x,y,z)
R = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z/R);
phi = atan(y./x);
end