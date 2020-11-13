function osn_pdf()
% PDF - probability distribudion function
t=tic;
global SPHERE_RADIUS D_THETA D_PHI ETTA D_T SHARP
N = 100;
c1_const = 0.0001;
SPHERE_RADIUS = 1.0;
D_THETA = c1_const*pi/180;
D_PHI = c1_const*pi/180;
D_T = .04*N^(-2);
ETTA = 1;
SHARP = 1;

[particles, duplicates] = GeneratePoints(N);
h = GenerateSphere(SPHERE_RADIUS,0,0,0);
h = VisualizePoints([particles duplicates],h);

for iter = 0:10000
    % Find potentials for each particle and duplicate
    maxVelocity = 0;
    sumPsi = 0;
    for k = 1:N
        % Current objects:
        particle = particles(k);
        duplicate = duplicates(k);
        
        % For the "other objects" choose all objects except the particle and
        % the duplicate number "k"
        otherIndices = ones(1,N);
        otherIndices(k) = 0;
        otherIndices = logical(otherIndices);
        otherObjects = [particles(otherIndices) duplicates(otherIndices)];
        
        %
        [particles(k).velocity, sumPsi1] = FindVelocity(particle, otherObjects);
        [duplicates(k).velocity, sumPsi2] = FindVelocity(duplicate, otherObjects);
        
        sumPsi = sumPsi + sumPsi1 + sumPsi2;        
        v2 = particles(k).velocity.theta.^2 + particles(k).velocity.phi.^2;
        maxVelocity = max([v2 maxVelocity]);
    end
   
    if maxVelocity < 3e-4
        D_T = .04*N^(-2)*(maxVelocity/3e-4)^(-1/2);
    end
    
    if maxVelocity < 3e-10
        break;
    end
    
    % Move each object according to its velocity
    for k = 1:N
       particles(k) = MoveObject(particles(k));
       duplicates(k) = MoveObject(duplicates(k));
    end
    
    if iter > 260 && rem(iter,20) == 0
        SHARP=SHARP+1;
    end
    
    if rem(iter,5) == 0
        h = VisualizePoints([particles duplicates], h);
        pause(0.01);
    end
    
    if rem(iter,1) == 0
        
        aTheta = [particles.theta duplicates.theta];
        aPhi = [particles.phi duplicates.phi];
        ins = 0;
        for k = 1:2*N
            ins = ins + IsInCone(aTheta(k), aPhi(k));
        end
        
        disp([...
            'iter: ' num2str(iter) ...
            ', maxVelocity: ' num2str(maxVelocity) ...
            ', sumPsi: ' num2str(sumPsi)...
            ', inside: ' num2str(ins)...
            ', outside: ' num2str(2*N - ins)]);
    end
    
    if rem(iter,40) == 0
        SaveData(particles);
    end
end
tt = toc(t);
disp(tt);
SaveData(particles);

end