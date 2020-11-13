function PDF_Check_PDFmatrix(aTheta, aPhi)
load('osn_pdf_output.mat');
oy = [0 1 0];
oz = [0 0 1];

sectorWidth = 1 / 180 * pi;
aSectors = sectorWidth:sectorWidth:pi/2;


v = [0 0 0];
for k = 1:length(aTheta)
    [v(1), v(2), v(3)] = Sph2Cart(1,aTheta(k),aPhi(k));
    
    alfa = acos(dot(oz, v)/norm(v));
    beta = acos(dot(oy, v)/norm(v));
        
    if alfa < pi/2 && beta < pi/2
        aAlfas = [aAlfas alfa];
    end
end

aAlfas = sort(aAlfas);
aBetas = sort(aBetas);



% Distribute angles by sectors
nSectors = length(aSectors);
aProbab.alfas = zeros(1,nSectors);
aProbab.betas = zeros(1,nSectors);
for sectorIdx = 1:nSectors
    aProbab.alfas(sectorIdx) = length(find(aAlfas < aSectors(sectorIdx))) / nAlfas;
    aProbab.betas(sectorIdx) = length(find(aBetas < aSectors(sectorIdx))) / nBetas;
end

% [exact.aTheta, exact.aI] = ExactSolution_PDFmatrix();

subplot(211);
plot(aSectors * 180 / pi,aProbab.alfas,'.-'); hold on;
plot(exact.aTheta * 180 / pi, exact.aI, 'r-');
grid;

subplot(212);
plot(aSectors(1:end-1) * 180 / pi,diff(smooth(aProbab,15)),'.-');
grid;
end
