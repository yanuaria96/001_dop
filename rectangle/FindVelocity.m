function [vx, vy, sumPsi] = FindVelocity(x_cur, y_cur, x_other, y_other)
% intput: (x_cur,y_cur) - current point
%         x_other, y_other - are x and y arrays for other points
% output: (vx,vy) - calculated velocity
%         sumPsi - summarized energy
global D_X D_Y ETTA

N = length(x_other);
sumPsi = 0;
vx = 0;
vy = 0;
for k = 1:N
    ox = x_other(k); % other x
    oy = y_other(k); % other y
    
    psiUnmoved = PSI(x_cur, y_cur, ox, oy);
    
    dPsidX = -(1/ETTA)*(1/D_X)*(PSI(x_cur+D_X, y_cur, ox, oy) - psiUnmoved);
    dPsidY = -(1/ETTA)*(1/D_Y)*(PSI(x_cur, y_cur+D_Y, ox, oy) - psiUnmoved);

    vx = vx + dPsidX;
    vy = vy + dPsidY;

    sumPsi = sumPsi + psiUnmoved;
end
PsiPenUnmoved = PsiPenalty(x_cur, y_cur);
dPsidX = -(1/ETTA)*(1/D_X)*(PsiPenalty(x_cur+D_X, y_cur) - PsiPenUnmoved);
dPsidY = -(1/ETTA)*(1/D_Y)*(PsiPenalty(x_cur, y_cur+D_Y) - PsiPenUnmoved);
vx = vx + dPsidX;
vy = vy + dPsidY;
sumPsi = sumPsi + PsiPenUnmoved*2;
end