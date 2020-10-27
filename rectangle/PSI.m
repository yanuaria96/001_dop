
function psi = PSI(x1, y1, x2, y2)
% input: (x1, y1) - point #1, (x2, y2) - point #2
% output: psi - energy

global W_2 H_2

if(x1 < W_2)
    f1 = 2;
else
    f1 = 1;
end
if(x2 < W_2)
    f2 = 2;
else
    f2 = 1;
end

rij = sqrt((x2-x1)^2 + (y2-y1)^2);

psi = 1/(rij^2)/sqrt(f1*f2);
%psi = 1/(rij^2)/(f1*f2)^(1/4);
%psi = -log(rij)/sqrt(f1*f2);
end