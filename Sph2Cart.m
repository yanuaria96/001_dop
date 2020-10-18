function [x,y,z] = Sph2Cart(R,theta,phi)
x=R.*sin(theta).*cos(phi);
y=R.*sin(theta).*sin(phi);
z=R.*cos(theta);
end