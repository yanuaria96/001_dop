function f=Mises_Fisher(k,myu, x)
%von Mises-Fisher distribution function
%k - cocentration parameter
%myu - vector on the unit sphere = mean direction
% x - vector on the unit sphere
C=k/(2*pi*(exp(k)-exp(-k)))*4*pi;
f= C*(exp(k*myu'*x)+exp(-k*myu'*x))/2;  %symmetrized Mises_Fisher
%f= C*exp(k*myu'*x); %classical Mises_Fisher
end