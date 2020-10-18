function [aTheta, I] = ExactSolution_PDFmatrix
%output I(theta), theta [0,pi]
myu = [0 0 1];
k = 2;

radius = 1;

nPhi = 600;
nTheta = 300;
dPhi = 2 * pi / nPhi;
dTheta = pi / 2/nTheta;

iSum = 0;
for n = 1:nTheta
    theta = n * dTheta;
    jSum = 0;
    for m = 1:nPhi
        phi = m * dPhi;
        [x(1),x(2),x(3)] = Sph2Cart(radius, theta, phi);
        %f(m) = Mises_Fisher(k, myu', x');
        %f(m) = Valee_Poussin(k,myu', x');
        f = PDF_Matrix(x');
        M = f*sin(theta);
        jSum = jSum+M;
    end 
    iSum = iSum + jSum;
    I(n) = iSum * dPhi * dTheta;
end
I = I/max(I);
plot([1:nTheta], I);

aTheta = [1:nTheta]*dTheta;
end
