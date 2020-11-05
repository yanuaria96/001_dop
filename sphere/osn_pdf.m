function osn_pdf()
% PDF - probability distribudion function
t=tic;
global SPHERE_RADIUS D_THETA D_PHI ETTA D_T
N = 20;
c1_const = 0.0001;
SPHERE_RADIUS = 1.0;
D_THETA = c1_const*pi/180;
D_PHI = c1_const*pi/180;
D_T = .04*N^(-2);
ETTA = 1;

[particles, duplicates] = GeneratePoints(N);
h = GenerateSphere(SPHERE_RADIUS,0,0,0);
h = VisualizePoints([particles duplicates],h);

for iter = 1:10000
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
    
    if rem(iter,1) == 0
        h = VisualizePoints([particles duplicates], h);
        pause(0.01);
    end
    
    if rem(iter,1) == 0
        disp([...
            'iter: ' num2str(iter) ...
            ', maxVelocity: ' num2str(maxVelocity) ...
            ', sumPsi: ' num2str(sumPsi)]);
    end
    
    if rem(iter,1) == 0
        SaveData(particles);
    end
end
tt = toc(t);
disp(tt);
SaveData(particles);

end