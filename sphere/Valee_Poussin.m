function f=Valee_Poussin(k,myu, x)
%von Mises-Fisher distribution function
%k - cocentration parameter
%myu - vector on the unit sphere = mean direction
% x - vector on the unit sphere
Omega = acos(dot(myu, x)/norm(x)/norm(myu));
Omega2 = pi - Omega;

f = (1/4)*beta(1.5,.5)/beta(1.5,k+.5) *((cos(Omega/2))^(2*k)+(cos(Omega2/2))^(2*k))/2;

end
