function [aTheta, I] = ExactSolution
myu = [0 0 1];
k = 2;

phi=0;

radius = 1;

nPoints = 200;
aTheta = linspace(0,90,nPoints) * pi / 180;



f = zeros(1,nPoints);
I = zeros(1,nPoints);
for m = 1:nPoints
    [x(1),x(2),x(3)] = Sph2Cart(radius, aTheta(m), phi);
    %f(m) = Mises_Fisher(k, myu', x');
    f(m) = Valee_Poussin(k,myu', x');
    %f(m) = PDF_Matrix(x');
    dI = f(1:m) .* sin(aTheta(1:m));
    I(m) = Integr(aTheta(1:m), dI);
end

I = I ./ max(I);
end

