function PDF_Check(aTheta, aPhi)
load('osn_pdf_output.mat');
oy = [0 0 1];
N = length(aTheta);
aAlfas = [];
v = [0 0 0];
for k = 1:N
    [v(1), v(2), v(3)] = Sph2Cart(1,aTheta(k),aPhi(k));
    alfa = acos(dot(oy, v)/norm(v));
    if alfa > pi/2
        continue;
    end
    aAlfas = [aAlfas alfa];
end
aAlfas = sort(aAlfas);
N2 = length(aAlfas);
sectorWidth = .5 / 180 * pi;
aSectors = sectorWidth:sectorWidth:pi/2;

% Distribute angles by sectors
nSectors = length(aSectors);
aProbab = zeros(1,nSectors);
sectorIdx = 1;
for k = 1:nSectors
    aProbab(sectorIdx) = length(find(aAlfas < aSectors(sectorIdx))) / N2;
    sectorIdx = sectorIdx + 1;
end


% 
[exact.aTheta, exact.aI] = ExactSolution();
%[exact.aTheta, exact.aI] = ExactSolution_PDFmatrix();

%subplot(211);
plot(aSectors * 180 / pi,aProbab,'.-'); hold on;
plot(exact.aTheta * 180 / pi, exact.aI, 'r-');
grid;

% subplot(212);
% plot(aSectors(1:end-1) * 180 / pi,diff(smooth(aProbab,15)),'.-');
% grid;
end
