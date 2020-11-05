
function [velocity, sumPsi] = FindVelocity(object, otherObjects)
global D_THETA D_PHI ETTA D_T 
N = length(otherObjects);
th1 = object.theta;
ph1 = object.phi;
velocity.theta = 0;
velocity.phi = 0;
sumPsi = 0;
for k = 1:N
    th2 = otherObjects(k).theta;
    ph2 = otherObjects(k).phi;
    
    psiUnmoved = PSI(th1,ph1,th2,ph2);
    
    dPsidTheta = - (D_T/ETTA).*(1/D_THETA).*(PSI(th1+D_THETA,ph1,th2,ph2) - psiUnmoved);
    velocity.theta = velocity.theta + dPsidTheta;
    
    dPsidPhi = - (D_T/ETTA).*(1/D_PHI).*(PSI(th1,ph1+D_PHI,th2,ph2)-psiUnmoved);    
    velocity.phi = velocity.phi + dPsidPhi;
    
    sumPsi = sumPsi + psiUnmoved;
end
end
