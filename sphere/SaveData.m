
function SaveData(points)
aTheta = [];
aPhi = [];
aX = [];
aY = [];
aZ = [];
for k = 1:length(points)
   aTheta(k) = points(k).theta;
   aPhi(k) = points(k).phi;
   [aX(k), aY(k), aZ(k)] = Sph2Cart(1, aTheta(k), aPhi(k));
end
dlmwrite('200_anis_x_y_z.txt', [aX' aY' aZ'],'delimiter','\t');
save 'osn_pdf_output';
end