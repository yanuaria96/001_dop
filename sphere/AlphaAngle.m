function alfa = AlphaAngle(theta, phi) 
% input: theta, phi - coordinates of the particle (1x1, 1x1)
% output: alfa (in radians)

[v(1), v(2), v(3)] = Sph2Cart(1, theta, phi);

% Find an angle between oy and vector v in 3D
oy = [0 1 0];
alfa = acos(dot(oy, v)/norm(v));
end