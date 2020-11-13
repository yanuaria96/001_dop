clc;
clear;

Pmax=2;
Pmin=1;
for i = 1:200
    beta(i) = 90*i/200;
    PDF(i) = 0.5*(Pmax+Pmin)+(Pmax-Pmin)/pi*atan(5*(beta(i)-30));

end
plot(beta, PDF);