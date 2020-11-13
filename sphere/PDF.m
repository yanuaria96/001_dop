function PDF = PDF(theta, phi)
% input: theta, phi - coordinates of the particle (1x1, 1x1) in radians
% output: pdf (1x1)
global SHARP

Pmax=2;
Pmin=1;
alfa = AlphaAngle(theta, phi);
beta = abs(alfa*180/pi -90);

% assert(beta >= 0);
% assert(beta <= 90);

PDF = 0.5*(Pmax+Pmin)+(Pmax-Pmin)/pi*atan(SHARP*(beta-30));

end